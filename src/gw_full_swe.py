#!/usr/bin/env python3
"""Full nonlinear pseudo-spectral shallow-water solver on an f-plane."""

from __future__ import annotations
import argparse
import math
from pathlib import Path
from typing import Any, Mapping
import numpy as np
from mod_bnd import make_sponge
from mod_ini import make_base_fplane, make_initial_ring_perturbation
from mod_spectral import (
    SpectralGrid,
    apply_spectral_filter,
    bilinear_sample,
    build_grid,
    compute_total_energy,
    load_config,
    resolve_output_directory,
    rk4_step,
    save_snapshot,
    spectral_dx,
    spectral_dy,
    write_power_timeseries,
    write_total_energy_row,
)

Array = np.ndarray
State = tuple[Array, Array, Array]


def apply_dealiasing_to_perturbations(
    h: Array,
    u: Array,
    v: Array,
    h_base: Array,
    u_base: Array,
    v_base: Array,
    grid: SpectralGrid,
    cutoff_fraction: float,) -> State:
    """Filter perturbations, then reconstruct the full SWE fields."""
    eta, u_perturbation, v_perturbation = apply_spectral_filter(
        h - h_base,
        u - u_base,
        v - v_base,
        grid,
        cutoff_fraction,
    )
    return (
        h_base + eta,
        u_base + u_perturbation,
        v_base + v_perturbation,
    )


def rhs_full_swe(
    h: Array,
    u: Array,
    v: Array,
    h_base: Array,
    u_base: Array,
    v_base: Array,
    sponge: Array,
    gravity: float,
    coriolis_parameter: float,
    grid: SpectralGrid,) -> State:
    """Evaluate the full nonlinear shallow-water equations."""
    h_x = spectral_dx(h, grid)
    h_y = spectral_dy(h, grid)
    u_x = spectral_dx(u, grid)
    u_y = spectral_dy(u, grid)
    v_x = spectral_dx(v, grid)
    v_y = spectral_dy(v, grid)

    h_t = -spectral_dx(h * u, grid) - spectral_dy(h * v, grid)
    u_t = -(u * u_x + v * u_y) + coriolis_parameter * v - gravity * h_x
    v_t = -(u * v_x + v * v_y) - coriolis_parameter * u - gravity * h_y

    h_t -= sponge * (h - h_base)
    u_t -= sponge * (u - u_base)
    v_t -= sponge * (v - v_base)
    return h_t, u_t, v_t


def compute_power_flux_full_swe(
    h: Array,
    u: Array,
    v: Array,
    h_base: Array,
    u_base: Array,
    v_base: Array,
    radius: float,
    density: float,
    gravity: float,
    mean_depth: float,
    circle_points: int,
    grid: SpectralGrid,) -> float:
    """Compute outward perturbation power through a diagnostic circle."""
    theta = np.linspace(0.0, 2.0 * np.pi, circle_points, endpoint=False)
    x_circle = radius * np.cos(theta)
    y_circle = radius * np.sin(theta)
    eta_circle = bilinear_sample(h - h_base, x_circle, y_circle, grid)
    u_circle = bilinear_sample(u - u_base, x_circle, y_circle, grid)
    v_circle = bilinear_sample(v - v_base, x_circle, y_circle, grid)
    radial_velocity = u_circle * np.cos(theta) + v_circle * np.sin(theta)
    radial_flux = density * gravity * mean_depth * eta_circle * radial_velocity
    return float(radius * np.trapezoid(radial_flux, theta))


def _full_swe_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    """Return the SWE-specific configuration with a clear missing-section error."""
    section = config.get("full_swe")
    if not isinstance(section, Mapping):
        raise ValueError("config.yaml must contain a 'full_swe' mapping.")
    return section


def run_full_swe(
    config: Mapping[str, Any],
    config_directory: str | Path,) -> None:
    """Build and integrate one configured full nonlinear SWE simulation."""
    grid = build_grid(config["grid"])
    physics = config["physics"]
    vortex = config["vortex"]
    initial = config["initial_condition"]
    sponge_config = config["sponge"]
    diagnostics = config["diagnostics"]
    numerics = config["numerics"]
    time = config["time"]
    full_swe = _full_swe_section(config)
    output = config["output"]

    gravity = float(physics["gravity_m_s2"])
    mean_depth = float(physics["mean_depth_m"])
    maximum_wind = float(vortex["maximum_wind_m_s"])
    vortex_radius = float(vortex["radius_km"])*1000.0
    coriolis_parameter = float(full_swe["coriolis_parameter_s_1"])
    density = float(diagnostics["air_density_kg_m3"])
    dt = float(time["time_step_s"])
    end_time = float(time["end_time_s"])
    snapshot_interval = float(output["snapshot_interval_s"])
    minimum_depth = float(full_swe["minimum_depth_m"])
    reference_radius_factor = float(full_swe["reference_radius_factor"])

    positive_values = {
        "physics.gravity_m_s2": gravity,
        "physics.mean_depth_m": mean_depth,
        "diagnostics.air_density_kg_m3": density,
        "time.time_step_s": dt,
        "time.end_time_s": end_time,
        "output.snapshot_interval_s": snapshot_interval,
        "full_swe.minimum_depth_m": minimum_depth,
    }
    for name, value in positive_values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive.")
    if not np.isfinite(coriolis_parameter):
        raise ValueError("full_swe.coriolis_parameter_s_1 must be finite.")
    if not np.isfinite(reference_radius_factor) or reference_radius_factor <= 1.0:
        raise ValueError("full_swe.reference_radius_factor must be greater than one.")

    nsteps = int(round(end_time / dt))
    output_stride = int(round(snapshot_interval / dt))
    if not math.isclose(nsteps * dt, end_time, abs_tol=1.0e-9):
        raise ValueError("time.end_time_s must be a multiple of its time.time_step_s.")
    if not math.isclose(output_stride * dt, snapshot_interval, abs_tol=1.0e-9):
        raise ValueError(
            "output.snapshot_interval_s must be a multiple of its time.time_step_s."
        )

    mode = str(initial["mode"]).strip().lower()
    if mode != "ring":
        raise ValueError("swe_spectral.py currently requires initial_condition.mode='ring'.")
    configured_ring_radius = float(initial["ring_radius_km"]) * 1000.0
    ring_radius = configured_ring_radius if configured_ring_radius > 0.0 else vortex_radius
    ring_width = float(initial["ring_width_km"]) * 1000.0
    dynamic_height = float(initial["amplitude_factor"]) * maximum_wind**2
    fixed_amplitude = float(initial["amplitude_magnetude_m"])

    corner_radius = 0.5 * np.hypot(grid.length_x, grid.length_y)
    reference_radius = reference_radius_factor * corner_radius
    u_base, v_base, h_base, _ = make_base_fplane(
        grid,
        maximum_wind,
        vortex_radius,
        gravity,
        mean_depth,
        coriolis_parameter,
        reference_radius,
        profile=str(initial["profile"]),
        smooth_kind=str(initial["smooth_kind"]),
        alpha=float(initial["alpha"]),
        nr=initial["nr"],
    )
    eta_initial, u_initial, v_initial = make_initial_ring_perturbation(
        grid,
        ring_radius,
        ring_width,
        dynamic_height,
        fixed_amplitude,
        gravity,
    )
    h = h_base + eta_initial
    u = u_base + u_initial
    v = v_base + v_initial
    sponge = make_sponge(
        grid,
        float(sponge_config["width_km"]) * 1000.0,
        float(sponge_config["damping_timescale_s"]),
    )

    diagnostic_radius = float(diagnostics["power_radius_km"]) * 1000.0
    maximum_circle_radius = 0.5 * min(grid.length_x, grid.length_y)
    if diagnostic_radius <= 0.0 or diagnostic_radius >= maximum_circle_radius:
        raise ValueError(
            "diagnostics.power_radius_km must be positive and smaller than "
            "half the shortest domain."
        )
    circle_points = int(diagnostics["power_circle_points"])
    if circle_points < 4:
        raise ValueError("diagnostics.power_circle_points must be at least four.")

    filter_enabled = bool(numerics["spectral_filter_enabled"])
    filter_fraction = float(numerics["spectral_filter_fraction"])
    if not 0.0 < filter_fraction <= 1.0:
        raise ValueError("numerics.spectral_filter_fraction must lie in (0, 1].")

    output_dir = resolve_output_directory(
        str(output["directory"]), Path(config_directory)
    )
    snapshot_dir = output_dir / "snapshots"
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    power_filename = str(output["outname"])+f"_{maximum_wind}_{vortex["radius_km"]}.csv"
    if not power_filename or Path(power_filename).name != power_filename:
        raise ValueError("output.outname must be a plain filename.")
    power_path = output_dir / power_filename
    energy_path = output_dir / "total_energy_timeseries.csv"

    rhs = lambda hh, uu, vv: rhs_full_swe(
        hh,
        uu,
        vv,
        h_base,
        u_base,
        v_base,
        sponge,
        gravity,
        coriolis_parameter,
        grid,
    )
    power_rows: list[tuple[float, float]] = []
    save_snapshot(snapshot_dir, 0, 0.0, u, v, h, grid)
    initial_power = compute_power_flux_full_swe(
        h,
        u,
        v,
        h_base,
        u_base,
        v_base,
        diagnostic_radius,
        density,
        gravity,
        mean_depth,
        circle_points,
        grid,
    )
    power_rows.append((0.0, initial_power))
    initial_energy = compute_total_energy(
        h - mean_depth,
        u,
        v,
        gravity,
        mean_depth,
    )
    write_total_energy_row(
        energy_path,
        0.0,
        initial_energy,
        write_header=True,
    )

    print(f"Grid: {grid.nx} x {grid.ny}, dx = {grid.dx / 1000.0:.1f} km")
    print(
        f"Domain: {grid.length_x / 1000.0:.0f} km x "
        f"{grid.length_y / 1000.0:.0f} km"
    )
    print(f"Coriolis parameter: {coriolis_parameter:.3e} s^-1")
    print(f"Reference radius: {reference_radius / 1000.0:.1f} km")
    print(f"Time step: {dt:.1f} s, total steps: {nsteps}")
    print(f"Snapshot output every {output_stride} steps")
    print(f"Power radius: {diagnostic_radius / 1000.0:.1f} km")
    print(
        f"step = 0, t = 0.000 h, "
        f"power flux = {initial_power:.6e} W, "
        f"total energy = {initial_energy:.12e} J"
    )

    for step in range(1, nsteps + 1):
        h, u, v = rk4_step(h, u, v, dt, rhs)
        if filter_enabled:
            h, u, v = apply_dealiasing_to_perturbations(
                h,
                u,
                v,
                h_base,
                u_base,
                v_base,
                grid,
                filter_fraction,
            )
        h = np.maximum(h, minimum_depth)
        time_s = step * dt

        if step % output_stride == 0 or step == nsteps:
            save_snapshot(snapshot_dir, step, time_s, u, v, h, grid)
            power = compute_power_flux_full_swe(
                h,
                u,
                v,
                h_base,
                u_base,
                v_base,
                diagnostic_radius,
                density,
                gravity,
                mean_depth,
                circle_points,
                grid,
            )
            power_rows.append((time_s, power))
            total_energy = compute_total_energy(
                h - mean_depth,
                u,
                v,
                gravity,
                mean_depth,
            )
            write_total_energy_row(energy_path, time_s, total_energy)
            print(
                f"step = {step}, t = {time_s / 3600.0:.3f} h, "
                f"power flux = {power:.6e} W, "
                f"total energy = {total_energy:.12e} J"
            )

    write_power_timeseries(power_path, power_rows)
    print(f"Saved snapshots to: {snapshot_dir}")
    print(f"Saved power time series to: {power_path}")
    print(f"Saved total energy time series to: {energy_path}")


def main() -> None:
    """Run the full SWE solver from the command line."""
    parser = argparse.ArgumentParser(
        description="Run the full nonlinear pseudo-spectral f-plane SWE model."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).with_name("config.yaml"),
        help="YAML configuration file (default: src/config.yaml)",
    )
    arguments = parser.parse_args()
    configuration_path = arguments.config.expanduser().resolve()
    configuration = load_config(configuration_path)
    run_full_swe(configuration, configuration_path.parent)


if __name__ == "__main__":
    main()
