"""PyTorch and torchrun backend for the pseudo-spectral wave model.

The solver has no trainable parameters, so torch.nn.parallel.DistributedDataParallel
is not applicable. Instead, torchrun ranks share independent batched FFT work and
sum the completed fields with torch.distributed collectives.
"""

from __future__ import annotations

import math
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import torch
import torch.distributed as dist

from mod_spectral import (
    RankineBaseState,
    SpectralGrid,
    build_grid,
    compute_power_flux,
    make_initial_perturbation,
    make_rankine_base,
    make_sponge,
    resolve_output_directory,
    save_snapshot,
    write_total_energy_row,
    write_power_timeseries,
)

Tensor = torch.Tensor
TensorState = tuple[Tensor, Tensor, Tensor]


@dataclass(frozen=True)
class DistributedContext:
    """torchrun rank and device information."""

    rank: int
    local_rank: int
    world_size: int
    device: torch.device
    distributed: bool

    @property
    def is_root(self) -> bool:
        return self.rank == 0


@dataclass(frozen=True)
class TorchGrid:
    """Tensor form of the spectral grid plus its NumPy output coordinates."""

    numpy: SpectralGrid
    KX: Tensor
    KY: Tensor
    K: Tensor


@dataclass(frozen=True)
class TorchBaseState:
    """Tensor form of the balanced Rankine-vortex state."""

    u: Tensor
    v: Tensor
    h: Tensor
    eta: Tensor
    u_x: Tensor
    u_y: Tensor
    v_x: Tensor
    v_y: Tensor


def initialize_distributed(
    device_option: str,
    require_distributed: bool,
) -> DistributedContext:
    """Initialize a torchrun process group and select one device per rank."""
    rank = int(os.environ.get("RANK", "0"))
    local_rank = int(os.environ.get("LOCAL_RANK", "0"))
    world_size = int(os.environ.get("WORLD_SIZE", "1"))
    distributed = world_size > 1
    if require_distributed and not distributed:
        raise RuntimeError(
            "--ddp requires torchrun with more than one process. For example: "
            "torchrun --standalone --nproc-per-node=2 ..."
        )

    if device_option == "auto":
        device_type = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device_type = device_option
    if device_type == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError(
                "CUDA was requested but no GPU is visible. Run inside a GPU "
                "allocation or use --device cpu."
            )
        if local_rank >= torch.cuda.device_count():
            raise RuntimeError(
                f"LOCAL_RANK={local_rank} exceeds the "
                f"{torch.cuda.device_count()} visible GPU(s)."
            )
        torch.cuda.set_device(local_rank)
        device = torch.device("cuda", local_rank)
    elif device_type == "cpu":
        device = torch.device("cpu")
    else:
        raise ValueError("device must be 'auto', 'cpu', or 'cuda'.")

    if distributed:
        backend = "nccl" if device.type == "cuda" else "gloo"
        dist.init_process_group(backend=backend, init_method="env://")
    return DistributedContext(rank, local_rank, world_size, device, distributed)


def torch_dtype(dtype_name: str) -> torch.dtype:
    """Translate a configured floating-point precision to a torch dtype."""
    try:
        return {"float32": torch.float32, "float64": torch.float64}[dtype_name]
    except KeyError as error:
        raise ValueError("dtype must be 'float32' or 'float64'.") from error


def _tensor(array: np.ndarray, context: DistributedContext, dtype: torch.dtype) -> Tensor:
    return torch.as_tensor(array, device=context.device, dtype=dtype)


def make_torch_grid(
    grid: SpectralGrid,
    context: DistributedContext,
    dtype: torch.dtype,
) -> TorchGrid:
    """Copy Fourier coordinates to the selected accelerator."""
    return TorchGrid(
        numpy=grid,
        KX=_tensor(grid.KX, context, dtype),
        KY=_tensor(grid.KY, context, dtype),
        K=_tensor(grid.K, context, dtype),
    )


def make_torch_base(
    base: RankineBaseState,
    context: DistributedContext,
    dtype: torch.dtype,
) -> TorchBaseState:
    """Copy a NumPy base state to the selected accelerator."""
    return TorchBaseState(
        *(
            _tensor(field, context, dtype)
            for field in (
                base.u,
                base.v,
                base.h,
                base.eta,
                base.u_x,
                base.u_y,
                base.v_x,
                base.v_y,
            )
        )
    )


def _all_reduce_sum(value: Tensor, context: DistributedContext) -> Tensor:
    """Sum fields for which exactly one rank computed each work item."""
    if context.distributed:
        target = torch.view_as_real(value) if value.is_complex() else value
        dist.all_reduce(target, op=dist.ReduceOp.SUM)
    return value


def distributed_derivatives(
    fields: tuple[Tensor, ...],
    wavenumbers: tuple[Tensor, ...],
    context: DistributedContext,
) -> tuple[Tensor, ...]:
    """Compute batched spectral derivatives, partitioned by work item."""
    count = len(fields)
    local_indices = list(range(context.rank, count, context.world_size))
    result = torch.zeros(
        (count, *fields[0].shape),
        device=context.device,
        dtype=fields[0].dtype,
    )
    if local_indices:
        field_batch = torch.stack([fields[index] for index in local_indices])
        wave_batch = torch.stack([wavenumbers[index] for index in local_indices])
        derivative_batch = torch.fft.ifft2(
            1j * wave_batch * torch.fft.fft2(field_batch)
        ).real
        result[local_indices] = derivative_batch
    _all_reduce_sum(result, context)
    return tuple(result.unbind(0))


def update_torch_base_state(
    initial_base: TorchBaseState,
    eta: Tensor,
    u: Tensor,
    v: Tensor,
    grid: TorchGrid,
    context: DistributedContext,
) -> TorchBaseState:
    """Build the distributed RHS background from initial plus perturbation."""
    du_dx, du_dy, dv_dx, dv_dy = distributed_derivatives(
        (u, u, v, v),
        (grid.KX, grid.KY, grid.KX, grid.KY),
        context,
    )
    return TorchBaseState(
        u=initial_base.u + u,
        v=initial_base.v + v,
        h=initial_base.h + eta,
        eta=initial_base.eta + eta,
        u_x=initial_base.u_x + du_dx,
        u_y=initial_base.u_y + du_dy,
        v_x=initial_base.v_x + dv_dx,
        v_y=initial_base.v_y + dv_dy,
    )


def distributed_fft2(
    fields: tuple[Tensor, ...],
    context: DistributedContext,
) -> tuple[Tensor, ...]:
    """Compute forward FFTs, partitioned by field across ranks."""
    count = len(fields)
    local_indices = list(range(context.rank, count, context.world_size))
    complex_dtype = (
        torch.complex64 if fields[0].dtype == torch.float32 else torch.complex128
    )
    result = torch.zeros(
        (count, *fields[0].shape),
        device=context.device,
        dtype=complex_dtype,
    )
    if local_indices:
        field_batch = torch.stack([fields[index] for index in local_indices])
        result[local_indices] = torch.fft.fft2(field_batch)
    _all_reduce_sum(result, context)
    return tuple(result.unbind(0))


def distributed_ifft2_real(
    fields: tuple[Tensor, ...],
    context: DistributedContext,
) -> TensorState:
    """Compute inverse FFTs, partitioned by field across ranks."""
    count = len(fields)
    local_indices = list(range(context.rank, count, context.world_size))
    real_dtype = (
        torch.float32 if fields[0].dtype == torch.complex64 else torch.float64
    )
    result = torch.zeros(
        (count, *fields[0].shape),
        device=context.device,
        dtype=real_dtype,
    )
    if local_indices:
        field_batch = torch.stack([fields[index] for index in local_indices])
        result[local_indices] = torch.fft.ifft2(field_batch).real
    _all_reduce_sum(result, context)
    fields_out = tuple(result.unbind(0))
    return fields_out[0], fields_out[1], fields_out[2]


def rhs_explicit_torch(
    eta: Tensor,
    u: Tensor,
    v: Tensor,
    base: TorchBaseState,
    sponge: Tensor,
    grid: TorchGrid,
    context: DistributedContext,
) -> TensorState:
    """Evaluate base-flow coupling and sponge terms with parallel FFT work."""
    derivatives = distributed_derivatives(
        (
            base.eta * u + eta * base.u,
            base.eta * v + eta * base.v,
            u,
            u,
            v,
            v,
        ),
        (grid.KX, grid.KY, grid.KX, grid.KY, grid.KX, grid.KY),
        context,
    )
    flux_x_dx, flux_y_dy, du_dx, du_dy, dv_dx, dv_dy = derivatives
    deta_dt = -flux_x_dx - flux_y_dy - sponge * eta
    du_dt = -(
        base.u * du_dx + base.v * du_dy + u * base.u_x + v * base.u_y
    ) - sponge * u
    dv_dt = -(
        base.u * dv_dx + base.v * dv_dy + u * base.v_x + v * base.v_y
    ) - sponge * v
    return deta_dt, du_dt, dv_dt


def rk4_step_torch(
    eta: Tensor,
    u: Tensor,
    v: Tensor,
    dt_step: float,
    rhs,
) -> TensorState:
    """Advance an explicit substep with fourth-order Runge--Kutta."""
    k1 = rhs(eta, u, v)
    k2 = rhs(
        eta + 0.5 * dt_step * k1[0],
        u + 0.5 * dt_step * k1[1],
        v + 0.5 * dt_step * k1[2],
    )
    k3 = rhs(
        eta + 0.5 * dt_step * k2[0],
        u + 0.5 * dt_step * k2[1],
        v + 0.5 * dt_step * k2[2],
    )
    k4 = rhs(
        eta + dt_step * k3[0],
        u + dt_step * k3[1],
        v + dt_step * k3[2],
    )
    return (
        eta + dt_step * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6,
        u + dt_step * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6,
        v + dt_step * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6,
    )


def exact_linear_gravity_step_torch(
    eta: Tensor,
    u: Tensor,
    v: Tensor,
    dt_step: float,
    grid: TorchGrid,
    gravity: float,
    mean_depth: float,
    context: DistributedContext,
) -> TensorState:
    """Advance the exact gravity-wave substep with distributed FFT batches."""
    E, U, V = distributed_fft2((eta, u, v), context)
    mask = grid.K > 0
    safe_K = torch.where(mask, grid.K, torch.ones_like(grid.K))
    longitudinal = (grid.KX * U + grid.KY * V) / safe_K
    transverse = (-grid.KY * U + grid.KX * V) / safe_K
    omega_dt = math.sqrt(gravity * mean_depth) * grid.K * dt_step
    coswt = torch.cos(omega_dt)
    sinwt = torch.sin(omega_dt)
    E_new = (
        E * coswt
        - 1j * math.sqrt(mean_depth / gravity) * longitudinal * sinwt
    )
    longitudinal_new = (
        longitudinal * coswt
        - 1j * math.sqrt(gravity / mean_depth) * E * sinwt
    )
    U_reconstructed = (
        grid.KX * longitudinal_new - grid.KY * transverse
    ) / safe_K
    V_reconstructed = (
        grid.KY * longitudinal_new + grid.KX * transverse
    ) / safe_K
    U_new = torch.where(mask, U_reconstructed, U)
    V_new = torch.where(mask, V_reconstructed, V)
    return distributed_ifft2_real((E_new, U_new, V_new), context)


def apply_spectral_filter_torch(
    eta: Tensor,
    u: Tensor,
    v: Tensor,
    grid: TorchGrid,
    cutoff_fraction: float,
    context: DistributedContext,
) -> TensorState:
    """Apply the configured cutoff, partitioned by field across ranks."""
    if not 0 < cutoff_fraction <= 1:
        raise ValueError("spectral_filter_fraction must lie in (0, 1].")
    kx_cut = cutoff_fraction * torch.max(torch.abs(grid.KX))
    ky_cut = cutoff_fraction * torch.max(torch.abs(grid.KY))
    keep = (torch.abs(grid.KX) <= kx_cut) & (torch.abs(grid.KY) <= ky_cut)
    fields = (eta, u, v)
    local_indices = list(range(context.rank, 3, context.world_size))
    result = torch.zeros(
        (3, *eta.shape), device=context.device, dtype=eta.dtype
    )
    if local_indices:
        batch = torch.stack([fields[index] for index in local_indices])
        result[local_indices] = torch.fft.ifft2(
            torch.fft.fft2(batch) * keep
        ).real
    _all_reduce_sum(result, context)
    filtered = tuple(result.unbind(0))
    return filtered[0], filtered[1], filtered[2]


def tensors_to_numpy(*fields: Tensor) -> tuple[np.ndarray, ...]:
    """Copy output fields from an accelerator to NumPy arrays."""
    return tuple(field.detach().cpu().numpy() for field in fields)


def compute_total_energy_torch(
    eta: Tensor,
    u: Tensor,
    v: Tensor,
    gravity: float,
    mean_depth: float,
) -> float:
    """Sum kinetic and gravitational energy on the active torch device."""
    depth = eta + mean_depth
    energy = (0.5 * u**2 + 0.5 * v**2 + 0.5 * gravity * depth) * depth
    return float(torch.sum(energy).item())


def _run_torch_model(
    config: Mapping[str, Any],
    config_directory: Path,
    context: DistributedContext,
    dtype: torch.dtype,
) -> None:
    """Execute the configured solver after distributed initialization."""
    grid_np = build_grid(config["grid"])
    physics = config["physics"]
    vortex = config["vortex"]
    initial = config["initial_condition"]
    time_config = config["time"]
    sponge_config = config["sponge"]
    diagnostics = config["diagnostics"]
    numerics = config["numerics"]
    output = config["output"]

    gravity = float(physics["gravity_m_s2"])
    mean_depth = float(physics["mean_depth_m"])
    base_update = physics.get("base_update", False)
    if not isinstance(base_update, bool):
        raise ValueError("physics.base_update must be true or false.")
    maximum_wind = float(vortex["maximum_wind_m_s"])
    vortex_radius = float(vortex["radius_km"]) * 1000
    dt = float(time_config["time_step_s"])
    end_time = float(time_config["end_time_s"])
    snapshot_interval = float(output["snapshot_interval_s"])
    nsteps = int(round(end_time / dt))
    output_stride = int(round(snapshot_interval / dt))
    if min(dt, end_time, snapshot_interval) <= 0:
        raise ValueError("Time step, end time, and output interval must be positive.")
    if not math.isclose(nsteps * dt, end_time, abs_tol=1e-12):
        raise ValueError("end_time_s must be an integer multiple of time_step_s.")
    if not math.isclose(output_stride * dt, snapshot_interval, abs_tol=1e-12):
        raise ValueError("snapshot_interval_s must be a multiple of time_step_s.")

    power_radius = float(diagnostics["power_radius_km"]) * 1000
    if power_radius >= 0.5 * min(grid_np.length_x, grid_np.length_y):
        raise ValueError(
            "power_radius_km must be smaller than half the shortest domain."
        )
    circle_points = int(diagnostics["power_circle_points"])
    if circle_points < 4:
        raise ValueError("power_circle_points must be at least 4.")

    base_np = make_rankine_base(
        grid_np, maximum_wind, vortex_radius, gravity, mean_depth
    )
    dynamic_height = float(initial["amplitude_factor"]) * maximum_wind**2
    eta_np, u_np, v_np = make_initial_perturbation(
        grid_np,
        vortex_radius,
        float(initial["ring_width_km"]) * 1000,
        dynamic_height,
        gravity,
    )
    sponge_np = make_sponge(
        grid_np,
        float(sponge_config["width_km"]) * 1000,
        float(sponge_config["damping_timescale_s"]),
    )
    grid = make_torch_grid(grid_np, context, dtype)
    initial_base = make_torch_base(base_np, context, dtype)
    base = initial_base
    eta = _tensor(eta_np, context, dtype)
    u = _tensor(u_np, context, dtype)
    v = _tensor(v_np, context, dtype)
    sponge = _tensor(sponge_np, context, dtype)
    if base_update:
        base = update_torch_base_state(
            initial_base, eta, u, v, grid, context
        )

    output_dir = resolve_output_directory(
        str(output["directory"]), config_directory
    )
    snapshot_dir = output_dir / "snapshots"
    energy_path = output_dir / "total_energy_timeseries.csv"
    density = float(diagnostics["air_density_kg_m3"])
    power_rows: list[tuple[float, float]] = []
    if context.is_root:
        snapshot_dir.mkdir(parents=True, exist_ok=True)
        eta_out, u_out, v_out = tensors_to_numpy(eta, u, v)
        save_snapshot(
            snapshot_dir,
            0,
            0,
            base_np.u + u_out,
            base_np.v + v_out,
            base_np.h + eta_out,
            grid_np,
        )
        power = compute_power_flux(
            eta_out,
            u_out,
            v_out,
            power_radius,
            density,
            gravity,
            circle_points,
            grid_np,
            mean_depth,
        )
        power_rows.append((0, power))
        total_energy = compute_total_energy_torch(
            base.eta, base.u, base.v, gravity, mean_depth
        )
        write_total_energy_row(
            energy_path, 0.0, total_energy, write_header=True
        )
        print(
            f"Torch backend: device={context.device}, dtype={dtype}, "
            f"ranks={context.world_size}"
        )
        print(
            f"Grid: {grid_np.nx} x {grid_np.ny}, "
            f"dx = {grid_np.dx / 1000:.1f} km"
        )
        print(f"Time step: {dt:.1f} s, total steps: {nsteps}")
        print(
            "Dynamic RHS base update: "
            f"{'enabled' if base_update else 'disabled'}"
        )
        print(f"Initial power flux: {power:.6e} W")
        print(f"step = 0, t = 0.000 h, total energy = {total_energy:.12e}")

    rhs = lambda ee, uu, vv: rhs_explicit_torch(
        ee, uu, vv, base, sponge, grid, context
    )
    filter_enabled = bool(numerics["spectral_filter_enabled"])
    filter_fraction = float(numerics["spectral_filter_fraction"])
    start_time = time.perf_counter()

    with torch.inference_mode():
        for step in range(1, nsteps + 1):
            eta, u, v = rk4_step_torch(eta, u, v, 0.5 * dt, rhs)
            eta, u, v = exact_linear_gravity_step_torch(
                eta, u, v, dt, grid, gravity, mean_depth, context
            )
            eta, u, v = rk4_step_torch(eta, u, v, 0.5 * dt, rhs)
            if filter_enabled:
                eta, u, v = apply_spectral_filter_torch(
                    eta, u, v, grid, filter_fraction, context
                )

            time_s = step * dt
            if base_update:
                base = update_torch_base_state(
                    initial_base, eta, u, v, grid, context
                )
            if context.is_root:
                total_energy = compute_total_energy_torch(
                    base.eta, base.u, base.v, gravity, mean_depth
                )

            if context.is_root and (
                step % output_stride == 0 or step == nsteps
            ):
                eta_out, u_out, v_out = tensors_to_numpy(eta, u, v)
                save_snapshot(
                    snapshot_dir,
                    step,
                    time_s,
                    base_np.u + u_out,
                    base_np.v + v_out,
                    base_np.h + eta_out,
                    grid_np,
                )
                power = compute_power_flux(
                    eta_out,
                    u_out,
                    v_out,
                    power_radius,
                    density,
                    gravity,
                    circle_points,
                    grid_np,
                    mean_depth,
                )
                power_rows.append((time_s, power))
                write_total_energy_row(energy_path, time_s, total_energy)
                print(
                    f"t = {time_s / 3600:.1f} h, "
                    f"power flux = {power:.6e} W, "
                    f"total energy = {total_energy:.12e} J"
                )

    if context.distributed:
        dist.barrier()
    if context.device.type == "cuda":
        torch.cuda.synchronize(context.device)
    if context.is_root:
        elapsed = time.perf_counter() - start_time
        power_path = output_dir / "power_flux_timeseries.csv"
        write_power_timeseries(power_path, power_rows)
        print(f"Elapsed integration time: {elapsed:.2f} s")
        print(f"Saved snapshots to: {snapshot_dir}")
        print(f"Saved power time series to: {power_path}")
        print(f"Saved total energy time series to: {energy_path}")


def run_torch_model(
    config: Mapping[str, Any],
    config_directory: str | Path,
    device: str = "auto",
    dtype_name: str = "float64",
    require_distributed: bool = False,
) -> None:
    """Run on PyTorch, optionally sharing FFT work across torchrun ranks."""
    context = initialize_distributed(device, require_distributed)
    try:
        _run_torch_model(
            config,
            Path(config_directory),
            context,
            torch_dtype(dtype_name),
        )
    finally:
        if dist.is_initialized():
            dist.destroy_process_group()
