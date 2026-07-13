"""Numerical routines for the pseudo-spectral TC gravity-wave model."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np

Array = np.ndarray
State = tuple[Array, Array, Array]
RhsFunction = Callable[[Array, Array, Array], State]


@dataclass(frozen=True)
class SpectralGrid:
    """Physical and Fourier-space coordinates for the periodic model grid."""

    length_x: float
    length_y: float
    dx: float
    dy: float
    nx: int
    ny: int
    x: Array
    y: Array
    X: Array
    Y: Array
    kx_1d: Array
    ky_1d: Array
    KX: Array
    KY: Array
    K: Array


@dataclass(frozen=True)
class RankineBaseState:
    """Balanced Rankine-vortex fields and analytic velocity gradients."""

    u: Array
    v: Array
    h: Array
    eta: Array
    u_x: Array
    u_y: Array
    v_x: Array
    v_y: Array


def load_config(config_path: str | Path) -> dict[str, Any]:
    """Load a YAML model configuration and check its top-level sections."""
    import yaml

    path = Path(config_path)
    with path.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration must be a YAML mapping: {path}")

    required = {
        "grid", "physics", "vortex", "initial_condition", "time",
        "sponge", "diagnostics", "numerics", "output",
    }
    missing = sorted(required.difference(config))
    if missing:
        raise ValueError(f"Missing configuration section(s): {', '.join(missing)}")
    return config


def build_grid(grid_config: Mapping[str, float]) -> SpectralGrid:
    """Create the regular physical grid and its Fourier wavenumbers."""
    length_x = float(grid_config["length_x_km"]) * 1000.0
    length_y = float(grid_config["length_y_km"]) * 1000.0
    dx = float(grid_config["spacing_x_km"]) * 1000.0
    dy = float(grid_config["spacing_y_km"]) * 1000.0
    if min(length_x, length_y, dx, dy) <= 0.0:
        raise ValueError("Grid lengths and spacings must be positive.")

    nx = int(round(length_x / dx))
    ny = int(round(length_y / dy))
    if not math.isclose(nx * dx, length_x, abs_tol=1.0e-9):
        raise ValueError("grid.length_x_km / grid.spacing_x_km must be an integer.")
    if not math.isclose(ny * dy, length_y, abs_tol=1.0e-9):
        raise ValueError("grid.length_y_km / grid.spacing_y_km must be an integer.")

    x = np.arange(nx) * dx - 0.5 * length_x
    y = np.arange(ny) * dy - 0.5 * length_y
    X, Y = np.meshgrid(x, y)
    kx_1d = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky_1d = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)
    KX, KY = np.meshgrid(kx_1d, ky_1d)
    K = np.sqrt(KX**2 + KY**2)
    return SpectralGrid(
        length_x, length_y, dx, dy, nx, ny, x, y, X, Y,
        kx_1d, ky_1d, KX, KY, K,
    )


def spectral_dx(field: Array, grid: SpectralGrid) -> Array:
    """Return the x derivative of a field using a Fourier transform."""
    transformed = np.fft.fft2(field)
    return np.fft.ifft2(1j * grid.KX * transformed).real


def spectral_dy(field: Array, grid: SpectralGrid) -> Array:
    """Return the y derivative of a field using a Fourier transform."""
    transformed = np.fft.fft2(field)
    return np.fft.ifft2(1j * grid.KY * transformed).real


def update_base_state(
    initial_base: RankineBaseState,
    eta: Array,
    u: Array,
    v: Array,
    grid: SpectralGrid,
) -> RankineBaseState:
    """Build the RHS background from the initial base plus the perturbation.

    The perturbation is retained. Velocity gradients are rebuilt from the
    analytic initial-base gradients plus spectral perturbation gradients.
    """
    return RankineBaseState(
        u=initial_base.u + u,
        v=initial_base.v + v,
        h=initial_base.h + eta,
        eta=initial_base.eta + eta,
        u_x=initial_base.u_x + spectral_dx(u, grid),
        u_y=initial_base.u_y + spectral_dy(u, grid),
        v_x=initial_base.v_x + spectral_dx(v, grid),
        v_y=initial_base.v_y + spectral_dy(v, grid),
    )


def make_rankine_base(
    grid: SpectralGrid,
    maximum_wind: float,
    radius: float,
    gravity: float,
    mean_depth: float,
) -> RankineBaseState:
    """Construct a balanced Rankine vortex in Cartesian coordinates."""
    if radius <= 0.0 or gravity <= 0.0 or mean_depth <= 0.0:
        raise ValueError("Vortex radius, gravity, and mean depth must be positive.")
    X, Y = grid.X, grid.Y
    r = np.sqrt(X**2 + Y**2)
    r_safe = np.maximum(r, 1.0e-12)
    u_theta = np.where(
        r <= radius,
        maximum_wind * r / radius,
        maximum_wind * radius / r_safe,
    )
    u = -u_theta * Y / r_safe
    v = u_theta * X / r_safe
    eta = np.where(
        r <= radius,
        -(maximum_wind**2 / gravity) * (1.0 - 0.5 * (r / radius) ** 2),
        -(maximum_wind**2 * radius**2) / (2.0 * gravity * r_safe**2),
    )
    h = mean_depth + eta

    u_x = np.zeros_like(X)
    u_y = np.zeros_like(X)
    v_x = np.zeros_like(X)
    v_y = np.zeros_like(X)
    inside = r <= radius
    outside = ~inside
    u_y[inside] = -maximum_wind / radius
    v_x[inside] = maximum_wind / radius

    scale = maximum_wind * radius
    r4 = r_safe**4
    u_x[outside] = 2.0 * scale * X[outside] * Y[outside] / r4[outside]
    u_y[outside] = -scale * (X[outside]**2 - Y[outside]**2) / r4[outside]
    v_x[outside] = scale * (Y[outside]**2 - X[outside]**2) / r4[outside]
    v_y[outside] = -2.0 * scale * X[outside] * Y[outside] / r4[outside]
    return RankineBaseState(u, v, h, eta, u_x, u_y, v_x, v_y)


def make_initial_perturbation(
    grid: SpectralGrid,
    radius: float,
    ring_width: float,
    dynamic_height_amplitude: float,
    gravity: float,
) -> State:
    """Create a Gaussian free-surface ring with zero initial velocity."""
    if ring_width <= 0.0:
        raise ValueError("Initial-condition ring width must be positive.")
    r = np.sqrt(grid.X**2 + grid.Y**2)
    eta_amplitude = dynamic_height_amplitude / gravity
    eta = eta_amplitude * np.exp(-((r - radius) ** 2) / (2.0 * ring_width**2))
    return eta, np.zeros_like(eta), np.zeros_like(eta)


def make_sponge(
    grid: SpectralGrid,
    sponge_width: float,
    damping_timescale: float,
) -> Array:
    """Create a quadratic damping coefficient near all four boundaries."""
    if sponge_width <= 0.0 or damping_timescale <= 0.0:
        raise ValueError("Sponge width and damping timescale must be positive.")
    if sponge_width >= 0.5 * min(grid.length_x, grid.length_y):
        raise ValueError("Sponge width must be smaller than half the domain.")
    dist_to_edge = np.minimum.reduce([
        grid.X + 0.5 * grid.length_x,
        0.5 * grid.length_x - grid.X,
        grid.Y + 0.5 * grid.length_y,
        0.5 * grid.length_y - grid.Y,
    ])
    ramp = np.clip((sponge_width - dist_to_edge) / sponge_width, 0.0, 1.0)
    return ramp**2 / damping_timescale


def rk4_step(
    eta: Array, u: Array, v: Array, dt_step: float, rhs: RhsFunction
) -> State:
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


def exact_linear_gravity_step(
    eta: Array,
    u: Array,
    v: Array,
    dt_step: float,
    grid: SpectralGrid,
    gravity: float,
    mean_depth: float,
) -> State:
    """Advance the linear gravity-wave system exactly in Fourier space."""
    E = np.fft.fft2(eta)
    U = np.fft.fft2(u)
    V = np.fft.fft2(v)
    longitudinal = np.zeros_like(E)
    transverse = np.zeros_like(E)
    mask = grid.K > 0.0
    longitudinal[mask] = (
        grid.KX[mask] * U[mask] + grid.KY[mask] * V[mask]
    ) / grid.K[mask]
    transverse[mask] = (
        -grid.KY[mask] * U[mask] + grid.KX[mask] * V[mask]
    ) / grid.K[mask]

    omega = math.sqrt(gravity * mean_depth) * grid.K
    coswt = np.cos(omega * dt_step)
    sinwt = np.sin(omega * dt_step)
    E_new = E.copy()
    longitudinal_new = longitudinal.copy()
    E_new[mask] = (
        E[mask] * coswt[mask]
        - 1j * math.sqrt(mean_depth / gravity)
        * longitudinal[mask] * sinwt[mask]
    )
    longitudinal_new[mask] = (
        longitudinal[mask] * coswt[mask]
        - 1j * math.sqrt(gravity / mean_depth) * E[mask] * sinwt[mask]
    )

    U_new = U.copy()
    V_new = V.copy()
    U_new[mask] = (
        grid.KX[mask] * longitudinal_new[mask]
        - grid.KY[mask] * transverse[mask]
    ) / grid.K[mask]
    V_new[mask] = (
        grid.KY[mask] * longitudinal_new[mask]
        + grid.KX[mask] * transverse[mask]
    ) / grid.K[mask]
    return (
        np.fft.ifft2(E_new).real,
        np.fft.ifft2(U_new).real,
        np.fft.ifft2(V_new).real,
    )


def rhs_explicit(
    eta: Array,
    u: Array,
    v: Array,
    base: RankineBaseState,
    sponge: Array,
    grid: SpectralGrid,
) -> State:
    """Evaluate base-flow coupling and sponge terms."""
    deta_dt = (
        -spectral_dx(base.eta * u + eta * base.u, grid)
        - spectral_dy(base.eta * v + eta * base.v, grid)
        - sponge * eta
    )
    du_dx = spectral_dx(u, grid)
    du_dy = spectral_dy(u, grid)
    dv_dx = spectral_dx(v, grid)
    dv_dy = spectral_dy(v, grid)
    du_dt = -(
        base.u * du_dx + base.v * du_dy + u * base.u_x + v * base.u_y
    ) - sponge * u
    dv_dt = -(
        base.u * dv_dx + base.v * dv_dy + u * base.v_x + v * base.v_y
    ) - sponge * v
    return deta_dt, du_dt, dv_dt


def apply_spectral_filter(
    eta: Array,
    u: Array,
    v: Array,
    grid: SpectralGrid,
    cutoff_fraction: float,
) -> State:
    """Apply a rectangular spectral cutoff to all perturbation fields."""
    if not 0.0 < cutoff_fraction <= 1.0:
        raise ValueError("spectral_filter_fraction must lie in (0, 1].")
    kx_cut = cutoff_fraction * np.max(np.abs(grid.kx_1d))
    ky_cut = cutoff_fraction * np.max(np.abs(grid.ky_1d))
    keep = (np.abs(grid.KX) <= kx_cut) & (np.abs(grid.KY) <= ky_cut)
    result = []
    for field in (eta, u, v):
        transformed = np.fft.fft2(field)
        transformed *= keep
        result.append(np.fft.ifft2(transformed).real)
    return result[0], result[1], result[2]


def bilinear_sample(
    field: Array, x_query: Array, y_query: Array, grid: SpectralGrid
) -> Array:
    """Bilinearly interpolate a regular-grid field inside the domain."""
    ix = (x_query - grid.x[0]) / grid.dx
    iy = (y_query - grid.y[0]) / grid.dy
    i0 = np.clip(np.floor(ix).astype(int), 0, field.shape[1] - 2)
    j0 = np.clip(np.floor(iy).astype(int), 0, field.shape[0] - 2)
    fx = ix - i0
    fy = iy - j0
    return (
        (1 - fx) * (1 - fy) * field[j0, i0]
        + fx * (1 - fy) * field[j0, i0 + 1]
        + (1 - fx) * fy * field[j0 + 1, i0]
        + fx * fy * field[j0 + 1, i0 + 1]
    )


def compute_power_flux(
    eta: Array,
    u: Array,
    v: Array,
    radius: float,
    density: float,
    gravity: float,
    circle_points: int,
    grid: SpectralGrid,
    H: float,
) -> float:
    """Compute signed outward perturbation power through a circle."""
    theta = np.linspace(0.0, 2.0 * np.pi, circle_points, endpoint=False)
    x_circle = radius * np.cos(theta)
    y_circle = radius * np.sin(theta)
    eta_circle = bilinear_sample(eta, x_circle, y_circle, grid)
    u_circle = bilinear_sample(u, x_circle, y_circle, grid)
    v_circle = bilinear_sample(v, x_circle, y_circle, grid)
    radial_velocity = u_circle * np.cos(theta) + v_circle * np.sin(theta)
    radial_flux = density * gravity * eta_circle * radial_velocity * H
    return float(radius * np.trapz(radial_flux, theta))


def compute_total_energy(
    eta: Array,
    u: Array,
    v: Array,
    gravity: float,
    mean_depth: float,
) -> float:
    """Sum kinetic and gravitational energy over all grid points."""
    depth = eta + mean_depth
    energy = (0.5 * u**2 + 0.5 * v**2 + 0.5 * gravity * depth) * depth
    return float(np.sum(energy, dtype=np.float64))


def save_snapshot(
    snapshot_dir: Path,
    step: int,
    time_s: float,
    u: Array,
    v: Array,
    h: Array,
    grid: SpectralGrid,
) -> Path:
    """Save one model snapshot as a compressed NumPy archive."""
    path = snapshot_dir / (
        f"snapshot_{step:06d}_t{int(round(time_s / 3600)):04d}h.npz"
    )
    np.savez_compressed(
        path,
        x=grid.x.astype(np.float32),
        y=grid.y.astype(np.float32),
        u=u.astype(np.float32),
        v=v.astype(np.float32),
        h=h.astype(np.float32),
        t=np.array([time_s], dtype=np.float64),
    )
    return path


def write_power_timeseries(
    path: Path, rows: list[tuple[float, float]]
) -> None:
    """Write power-flux diagnostic rows to CSV."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", "power_W"])
        writer.writerows(
            (f"{time_s:.6f}", f"{power:.12e}") for time_s, power in rows
        )


def write_total_energy_row(
    path: Path,
    time_s: float,
    total_energy: float,
    write_header: bool = False,
) -> None:
    """Write one energy row immediately so progress survives interrupted runs."""
    mode = "w" if write_header else "a"
    with path.open(mode, newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        if write_header:
            writer.writerow(["time_s", "total_energy"])
        writer.writerow([f"{time_s:.6f}", f"{total_energy:.12e}"])


def resolve_output_directory(
    output_value: str, config_directory: Path
) -> Path:
    """Resolve an output path relative to the YAML file."""
    output_dir = Path(output_value).expanduser()
    if not output_dir.is_absolute():
        output_dir = config_directory / output_dir
    return output_dir.resolve()


def run_model(config: Mapping[str, Any], config_directory: str | Path) -> None:
    """Build, integrate, diagnose, and save one model simulation."""
    grid = build_grid(config["grid"])
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
    vortex_radius_km = float(vortex["radius_km"])
    vortex_radius = vortex_radius_km * 1000.0
    dt = float(time_config["time_step_s"])
    end_time = float(time_config["end_time_s"])
    snapshot_interval = float(output["snapshot_interval_s"])
    if min(dt, end_time, snapshot_interval) <= 0.0:
        raise ValueError("Time step, end time, and output interval must be positive.")

    nsteps = int(round(end_time / dt))
    output_stride = int(round(snapshot_interval / dt))
    if not math.isclose(nsteps * dt, end_time, abs_tol=1.0e-12):
        raise ValueError("end_time_s must be an integer multiple of time_step_s.")
    if not math.isclose(
        output_stride * dt, snapshot_interval, abs_tol=1.0e-12
    ):
        raise ValueError("snapshot_interval_s must be a multiple of time_step_s.")

    power_radius = float(diagnostics["power_radius_km"]) * 1000.0
    if power_radius >= 0.5 * min(grid.length_x, grid.length_y):
        raise ValueError(
            "power_radius_km must be smaller than half the shortest domain."
        )

    initial_base = make_rankine_base(
        grid, maximum_wind, vortex_radius, gravity, mean_depth
    )
    base = initial_base
    dynamic_height = float(initial["amplitude_factor"]) * maximum_wind**2
    eta, u, v = make_initial_perturbation(
        grid,
        vortex_radius,
        float(initial["ring_width_km"]) * 1000.0,
        dynamic_height,
        gravity,
    )
    sponge = make_sponge(
        grid,
        float(sponge_config["width_km"]) * 1000.0,
        float(sponge_config["damping_timescale_s"]),
    )
    if base_update:
        base = update_base_state(initial_base, eta, u, v, grid)

    output_dir = resolve_output_directory(
        str(output["directory"]), Path(config_directory)
    )
    snapshot_dir = output_dir / "snapshots"
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    energy_path = output_dir / "total_energy_timeseries.csv"
    density = float(diagnostics["air_density_kg_m3"])
    circle_points = int(diagnostics["power_circle_points"])
    if circle_points < 4:
        raise ValueError("power_circle_points must be at least 4.")

    save_snapshot(
        snapshot_dir,
        0,
        0.0,
        initial_base.u + u,
        initial_base.v + v,
        initial_base.h + eta,
        grid,
    )
    power = compute_power_flux(
        eta, u, v, power_radius, density, gravity, circle_points, grid, mean_depth
    )
    power_rows = [(0.0, power)]
    total_energy = compute_total_energy(
        base.eta, base.u, base.v, gravity, mean_depth
    )
    write_total_energy_row(
        energy_path, 0.0, total_energy, write_header=True
    )
    print(f"Grid: {grid.nx} x {grid.ny}, dx = {grid.dx / 1000:.1f} km")
    print(
        f"Domain: {grid.length_x / 1000:.0f} km x "
        f"{grid.length_y / 1000:.0f} km"
    )
    print(f"Time step: {dt:.1f} s, total steps: {nsteps}")
    print(f"Snapshot output every {output_stride} steps")
    print(f"Power radius: {power_radius / 1000:.1f} km")
    print(f"Dynamic RHS base update: {'enabled' if base_update else 'disabled'}")
    print(f"Initial power flux: {power:.6e} W")
    print(f"step = 0, t = 0.000 h, total energy = {total_energy:.12e}")

    rhs: RhsFunction = lambda ee, uu, vv: rhs_explicit(
        ee, uu, vv, base, sponge, grid
    )
    filter_enabled = bool(numerics["spectral_filter_enabled"])
    filter_fraction = float(numerics["spectral_filter_fraction"])

    for step in range(1, nsteps + 1):
        eta, u, v = rk4_step(eta, u, v, 0.5 * dt, rhs)
        eta, u, v = exact_linear_gravity_step(
            eta, u, v, dt, grid, gravity, mean_depth
        )
        eta, u, v = rk4_step(eta, u, v, 0.5 * dt, rhs)
        if filter_enabled:
            eta, u, v = apply_spectral_filter(
                eta, u, v, grid, filter_fraction
            )

        time_s = step * dt
        if base_update:
            base = update_base_state(initial_base, eta, u, v, grid)
        total_energy = compute_total_energy(
            base.eta, base.u, base.v, gravity, mean_depth
        )

        if step % output_stride == 0 or step == nsteps:
            save_snapshot(
                snapshot_dir,
                step,
                time_s,
                initial_base.u + u,
                initial_base.v + v,
                initial_base.h + eta,
                grid,
            )
            power = compute_power_flux(
                eta, u, v, power_radius, density, gravity, circle_points, grid, mean_depth
            )
            power_rows.append((time_s, power))
            write_total_energy_row(energy_path, time_s, total_energy)
            print(
                f"t = {time_s / 3600:.3f} h, "
                f"power flux = {power:.6e} W, "
                f"total energy = {total_energy:.12e} J"
            )

    power_path = output_dir / (
        f"power_flux_timeseries_{maximum_wind}_{vortex_radius_km}.csv"
    )
    write_power_timeseries(power_path, power_rows)
    print(f"Saved snapshots to: {snapshot_dir}")
    print(f"Saved power time series to: {power_path}")
    print(f"Saved total energy time series to: {energy_path}")
