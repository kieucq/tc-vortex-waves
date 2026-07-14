"""Boundary damping routines for the pseudo-spectral model."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from mod_spectral import SpectralGrid

Array = np.ndarray


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
