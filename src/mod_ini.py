"""Initial-condition and balanced-vortex construction routines."""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterable

import numpy as np

if TYPE_CHECKING:
    from mod_spectral import RankineBaseState, SpectralGrid

State = tuple[np.ndarray, np.ndarray, np.ndarray]


def make_initial_ring_perturbation(
    grid: SpectralGrid,
    radius: float,
    ring_width: float,
    dynamic_height_amplitude: float,
    pertubation_amplitude: float,
    gravity: float) -> State:
    """Create a Gaussian free-surface ring with zero initial velocity."""
    if ring_width <= 0.0:
        raise ValueError("Initial-condition ring width must be positive.")
    r = np.sqrt(grid.X**2 + grid.Y**2)
    eta_amplitude = dynamic_height_amplitude / gravity
    if pertubation_amplitude > 1e-18:
         eta = pertubation_amplitude*np.exp(-((r - radius) ** 2) / (2.0 * ring_width**2))
    else:
         eta = eta_amplitude * np.exp(-((r - radius) ** 2) / (2.0 * ring_width**2))
    return eta, np.zeros_like(eta), np.zeros_like(eta)


def make_initial_vht_perturbation(
    model_time: float,
    perturbation_magnitude: float,
    perturbation_width_scale: float,
    locations: Iterable[tuple[float, float]],
    trigger_times: Iterable[float],
    grid: SpectralGrid,
) -> State:
    """Superpose geostrophically balanced Gaussian perturbations.

    ``locations`` contains the ``(x, y)`` center of each perturbation in the
    same coordinate units as the grid. ``perturbation_magnitude`` is the peak
    free-surface displacement, and ``perturbation_width_scale`` is the Gaussian
    standard deviation in the same length units as the grid. A location is
    activated only when its corresponding trigger time matches ``model_time``.
    """
    eta_amplitude = float(perturbation_magnitude)
    width_scale = float(perturbation_width_scale)
    if not np.isfinite(eta_amplitude):
        raise ValueError("Perturbation magnitude must be finite.")
    if not np.isfinite(width_scale) or width_scale <= 0.0:
        raise ValueError("Perturbation width scale must be finite and positive.")
    current_time = float(model_time)
    if not np.isfinite(current_time):
        raise ValueError("Model time must be finite.")
    spacing_x = float(grid.dx)
    spacing_y = float(grid.dy)
    if not np.isfinite(spacing_x) or spacing_x <= 0.0:
        raise ValueError("Grid spacing grid.dx must be finite and positive.")
    if not np.isfinite(spacing_y) or spacing_y <= 0.0:
        raise ValueError("Grid spacing grid.dy must be finite and positive.")

    centers = np.asarray(list(locations), dtype=float)
    if centers.size == 0:
        centers = np.empty((0, 2), dtype=float)
    elif centers.ndim == 1 and centers.shape == (2,):
        centers = centers.reshape(1, 2)
    if centers.ndim != 2 or centers.shape[1] != 2:
        raise ValueError("Locations must be an iterable of (x, y) pairs.")
    if not np.all(np.isfinite(centers)):
        raise ValueError("All perturbation locations must be finite.")
    times = np.asarray(list(trigger_times), dtype=float)
    if times.ndim != 1:
        raise ValueError("Trigger times must be a one-dimensional iterable.")
    if times.size != centers.shape[0]:
        raise ValueError("Locations and trigger times must have the same length.")
    if not np.all(np.isfinite(times)):
        raise ValueError("All trigger times must be finite.")
    if grid.X.ndim != 2 or min(grid.X.shape) < 2:
        raise ValueError("The grid must contain at least two points on each axis.")

    total_eta = np.zeros_like(grid.X, dtype=float)
    width_squared = width_scale**2

    active = np.isclose(times, current_time, rtol=0.0, atol=1.0e-9)
    for center_x, center_y in centers[active]:
        delta_x = grid.X - center_x
        delta_y = grid.Y - center_y
        radius = np.hypot(delta_x, delta_y)
        eta = eta_amplitude * np.exp(-(radius**2) / (2.0 * width_squared))
        total_eta += eta

    edge_order = 2 if min(total_eta.shape) >= 3 else 1
    eta_y, eta_x = np.gradient(
        total_eta,
        spacing_y,
        spacing_x,
        edge_order=edge_order,
    )
    gravity = 9.81
    coriolis_parameter = 1.0e-5
    geostrophic_factor = gravity / coriolis_parameter
    total_u = -geostrophic_factor * eta_y
    total_v = geostrophic_factor * eta_x

    return total_eta, total_u, total_v


def make_rankine_base(
    grid: SpectralGrid,
    maximum_wind: float,
    radius: float,
    gravity: float,
    mean_depth: float) -> RankineBaseState:
    """Construct a balanced Rankine vortex in Cartesian coordinates."""
    # Imported here to avoid a module-level cycle: mod_spectral imports this
    # initializer while also defining the shared base-state container.
    from mod_spectral import RankineBaseState

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
