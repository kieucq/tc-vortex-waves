#!/usr/bin/env python3
"""
Pseudo-spectral linear shallow-water solver about a balanced Rankine vortex.

Model:
    - Linearized SWEs about a steady Rankine vortex base state
    - Fourier pseudo-spectral derivatives on a periodic square grid
    - Exact gravity-wave substep in Fourier space
    - Explicit base-flow coupling substep with RK4
    - Sponge layer near boundaries to reduce reflection

Outputs:
    - Hourly snapshots of total (u, v, h) on the regular grid
    - Outward perturbation power flux through a circle at R_OUT

Notes:
    - This solves the perturbation equations about a fixed balanced Rankine vortex.
    - A true nonreflecting boundary condition is approximated by a sponge layer.
    - With L = 3000 km, the largest usable radius is < 1500 km.
"""

from __future__ import annotations

import os
import csv
import math
import numpy as np


# -----------------------------
# Physical and numerical inputs
# -----------------------------
Lx_km = 3000.0
Ly_km = 3000.0
dx_km = 10.0
dy_km = 10.0

Umax = 50.0          # m/s, max tangential wind at r = R
R_km = 50.0          # km
A_dyn = 0.1*Umax**2  # m^2/s^2, as given by the user
H = 5000.0           # m
g = 9.81             # m/s^2
dt = 90.0            # s
t_end = 3 * 24 * 3600.0  # 3 days

rho = 1000.0         # kg/m^3, water density for power diagnostic

# Default power-flux radius inside the domain.
# 2500 km would not fit in a 3000 km square; use a larger domain for that.
R_OUT_km = 1350.0

# Gaussian ring width
sigma_km = 20.0      # km          

# Sponge layer
sponge_width_km = 250.0
sponge_tau_s = 3600.0  # damping timescale near the boundary
sponge_max = 1.0 / sponge_tau_s

# Output
outdir = "outputs"
snapshot_dir = os.path.join(outdir, "snapshots")
os.makedirs(snapshot_dir, exist_ok=True)

# -----------------------------
# Grid
# -----------------------------
Lx = Lx_km * 1000.0
Ly = Ly_km * 1000.0
dx = dx_km * 1000.0
dy = dy_km * 1000.0

nx = int(round(Lx / dx))
ny = int(round(Ly / dy))

if abs(nx * dx - Lx) > 1e-9 or abs(ny * dy - Ly) > 1e-9:
    raise ValueError("Lx/dx and Ly/dy must be integers for this setup.")

x = np.arange(nx) * dx - 0.5 * Lx
y = np.arange(ny) * dy - 0.5 * Ly
X, Y = np.meshgrid(x, y)

x0 = x[0]
y0 = y[0]

# Fourier wavenumbers
kx_1d = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
ky_1d = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)
KX, KY = np.meshgrid(kx_1d, ky_1d)
K2 = KX**2 + KY**2
K = np.sqrt(K2)
c = math.sqrt(g * H)

# -----------------------------
# Helper functions
# -----------------------------
def spectral_dx(f: np.ndarray) -> np.ndarray:
    """Spectral x-derivative using FFT."""
    F = np.fft.fft2(f)
    return np.fft.ifft2(1j * KX * F).real


def spectral_dy(f: np.ndarray) -> np.ndarray:
    """Spectral y-derivative using FFT."""
    F = np.fft.fft2(f)
    return np.fft.ifft2(1j * KY * F).real


def make_rankine_base(X: np.ndarray, Y: np.ndarray, Umax: float, R: float, g: float, H: float):
    """
    Balanced Rankine vortex in Cartesian coordinates.

    Inside r <= R:
        u_theta = Umax * r / R
    Outside r > R:
        u_theta = Umax * R / r

    Cartesian velocity:
        u = -u_theta * y / r
        v =  u_theta * x / r

    Balanced free-surface displacement:
        eta0(r) = -Umax^2/g * (1 - r^2/(2R^2)),    r <= R
                = -Umax^2 R^2 / (2 g r^2),         r > R
    """
    eps = 1e-12
    r = np.sqrt(X**2 + Y**2)
    r_safe = np.maximum(r, eps)

    u_theta = np.where(r <= R, Umax * r / R, Umax * R / r_safe)

    u0 = -u_theta * Y / r_safe
    v0 =  u_theta * X / r_safe

    eta0 = np.where(
        r <= R,
        -(Umax**2 / g) * (1.0 - 0.5 * (r / R)**2),
        -(Umax**2 * R**2) / (2.0 * g * r_safe**2),
    )
    h0 = H + eta0

    # Analytic derivatives for the Cartesian Rankine vortex
    u0_x = np.zeros_like(X)
    u0_y = np.zeros_like(X)
    v0_x = np.zeros_like(X)
    v0_y = np.zeros_like(X)

    inside = r <= R
    outside = ~inside

    # Inside: solid-body rotation
    u0_x[inside] = 0.0
    u0_y[inside] = -Umax / R
    v0_x[inside] = Umax / R
    v0_y[inside] = 0.0

    # Outside: u = -U R y / r^2, v = U R x / r^2
    s = Umax * R
    r2 = r_safe**2
    r4 = r2**2

    u0_x[outside] = 2.0 * s * X[outside] * Y[outside] / r4[outside]
    u0_y[outside] = -s * (X[outside]**2 - Y[outside]**2) / r4[outside]
    v0_x[outside] = s * (Y[outside]**2 - X[outside]**2) / r4[outside]
    v0_y[outside] = -2.0 * s * X[outside] * Y[outside] / r4[outside]

    return u0, v0, h0, eta0, u0_x, u0_y, v0_x, v0_y


def make_sponge(X: np.ndarray, Y: np.ndarray, Lx: float, Ly: float, sponge_width: float, sigma_max: float) -> np.ndarray:
    """
    Quadratic sponge damping coefficient, zero in the interior and increasing
    smoothly toward the boundary.
    """
    x_min = -0.5 * Lx
    x_max =  0.5 * Lx
    y_min = -0.5 * Ly
    y_max =  0.5 * Ly

    dist_to_edge = np.minimum.reduce([
        X - x_min,
        x_max - X,
        Y - y_min,
        y_max - Y,
    ])

    ramp = np.clip((sponge_width - dist_to_edge) / sponge_width, 0.0, 1.0)
    return sigma_max * ramp**2


def rk4_step(eta: np.ndarray, u: np.ndarray, v: np.ndarray, dt_step: float,
             rhs_func) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Classic RK4 for the explicit coupling substep."""
    k1 = rhs_func(eta, u, v)
    k2 = rhs_func(eta + 0.5 * dt_step * k1[0],
                  u   + 0.5 * dt_step * k1[1],
                  v   + 0.5 * dt_step * k1[2])
    k3 = rhs_func(eta + 0.5 * dt_step * k2[0],
                  u   + 0.5 * dt_step * k2[1],
                  v   + 0.5 * dt_step * k2[2])
    k4 = rhs_func(eta + dt_step * k3[0],
                  u   + dt_step * k3[1],
                  v   + dt_step * k3[2])

    eta_new = eta + (dt_step / 6.0) * (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0])
    u_new   = u   + (dt_step / 6.0) * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])
    v_new   = v   + (dt_step / 6.0) * (k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2])

    return eta_new, u_new, v_new


def exact_linear_gravity_step(eta: np.ndarray, u: np.ndarray, v: np.ndarray, dt_step: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Exact gravity-wave substep in Fourier space for the constant-depth linear system:
        eta_t = -H (u_x + v_y)
        u_t   = -g eta_x
        v_t   = -g eta_y

    This stabilizes the fast gravity-wave part for dt = 90 s.
    """
    E = np.fft.fft2(eta)
    U = np.fft.fft2(u)
    V = np.fft.fft2(v)

    # Longitudinal / transverse decomposition in Fourier space
    up = np.zeros_like(E)
    ut = np.zeros_like(E)

    mask = K > 0.0
    up[mask] = (KX[mask] * U[mask] + KY[mask] * V[mask]) / K[mask]
    ut[mask] = (-KY[mask] * U[mask] + KX[mask] * V[mask]) / K[mask]

    omega = c * K
    coswt = np.cos(omega * dt_step)
    sinwt = np.sin(omega * dt_step)

    E_new = E.copy()
    up_new = up.copy()
    ut_new = ut.copy()

    # For K > 0, exact oscillator update.
    # E_new = E cos(wt) - i * sqrt(H/g) * up sin(wt)
    # up_new = up cos(wt) - i * sqrt(g/H) * E sin(wt)
    E_new[mask] = E[mask] * coswt[mask] - 1j * math.sqrt(H / g) * up[mask] * sinwt[mask]
    up_new[mask] = up[mask] * coswt[mask] - 1j * math.sqrt(g / H) * E[mask] * sinwt[mask]
    # transverse component unchanged
    ut_new[mask] = ut[mask]

    # Reconstruct Fourier-space U, V
    U_new = U.copy()
    V_new = V.copy()

    U_new[mask] = (KX[mask] * up_new[mask] - KY[mask] * ut_new[mask]) / K[mask]
    V_new[mask] = (KY[mask] * up_new[mask] + KX[mask] * ut_new[mask]) / K[mask]

    # k = 0 mode remains unchanged
    eta_new = np.fft.ifft2(E_new).real
    u_new = np.fft.ifft2(U_new).real
    v_new = np.fft.ifft2(V_new).real

    return eta_new, u_new, v_new


def rhs_explicit(eta: np.ndarray, u: np.ndarray, v: np.ndarray,
                 u0: np.ndarray, v0: np.ndarray, eta0: np.ndarray,
                 u0_x: np.ndarray, u0_y: np.ndarray, v0_x: np.ndarray, v0_y: np.ndarray,
                 sponge: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Explicit linearized coupling terms about the balanced Rankine vortex.

    Splitting:
        - exact linear gravity-wave part handled in Fourier space
        - remaining base-flow coupling handled explicitly

    Continuity explicit part:
        eta_t = -d_x(eta0 * u + eta * u0) - d_y(eta0 * v + eta * v0) - sponge * eta

    Momentum explicit part:
        u_t = -(u0 u_x + v0 u_y + u u0_x + v u0_y) - sponge * u
        v_t = -(u0 v_x + v0 v_y + u v0_x + v v0_y) - sponge * v
    """
    # continuity explicit term
    flux_x = eta0 * u + eta * u0
    flux_y = eta0 * v + eta * v0
    deta_dt = -spectral_dx(flux_x) - spectral_dy(flux_y) - sponge * eta

    # velocity derivatives
    du_dx = spectral_dx(u)
    du_dy = spectral_dy(u)
    dv_dx = spectral_dx(v)
    dv_dy = spectral_dy(v)

    du_dt = -(u0 * du_dx + v0 * du_dy + u * u0_x + v * u0_y) - sponge * u
    dv_dt = -(u0 * dv_dx + v0 * dv_dy + u * v0_x + v * v0_y) - sponge * v

    return deta_dt, du_dt, dv_dt


def bilinear_sample(field: np.ndarray, xq: np.ndarray, yq: np.ndarray, x0: float, y0: float, dx: float, dy: float) -> np.ndarray:
    """
    Bilinear interpolation on the regular grid for points known to lie inside the domain.
    field shape: (ny, nx)
    """
    ix = (xq - x0) / dx
    iy = (yq - y0) / dy

    i0 = np.floor(ix).astype(int)
    j0 = np.floor(iy).astype(int)

    # Clip so that i0+1 and j0+1 remain valid
    i0 = np.clip(i0, 0, field.shape[1] - 2)
    j0 = np.clip(j0, 0, field.shape[0] - 2)

    fx = ix - i0
    fy = iy - j0

    f00 = field[j0,     i0    ]
    f10 = field[j0,     i0 + 1]
    f01 = field[j0 + 1, i0    ]
    f11 = field[j0 + 1, i0 + 1]

    return ((1.0 - fx) * (1.0 - fy) * f00 +
            fx * (1.0 - fy) * f10 +
            (1.0 - fx) * fy * f01 +
            fx * fy * f11)


def compute_power_flux(eta_p: np.ndarray, u_p: np.ndarray, v_p: np.ndarray,
                       Rout: float, rho: float, g: float,
                       ntheta: int = 720) -> float:
    """
    Signed outward perturbation power flux across a circle of radius Rout:
        P = ∮ rho * g * eta' * u_r' * R dtheta
    """
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    xs = Rout * np.cos(theta)
    ys = Rout * np.sin(theta)

    eta_c = bilinear_sample(eta_p, xs, ys, x0, y0, dx, dy)
    u_c = bilinear_sample(u_p, xs, ys, x0, y0, dx, dy)
    v_c = bilinear_sample(v_p, xs, ys, x0, y0, dx, dy)

    ur_c = u_c * np.cos(theta) + v_c * np.sin(theta)
    fr = rho * g * eta_c * ur_c  # W/m^2
    P = Rout * np.trapz(fr, theta)  # Watts per unit depth (actually total across the circle)
    return float(P)


def save_snapshot(step: int, t: float, u_total: np.ndarray, v_total: np.ndarray, h_total: np.ndarray) -> None:
    """
    Save one hourly snapshot as a compressed NPZ file.
    """
    fname = os.path.join(snapshot_dir, f"snapshot_{step:06d}_t{int(round(t/3600)):04d}h.npz")
    np.savez_compressed(
        fname,
        x=x.astype(np.float32),
        y=y.astype(np.float32),
        u=u_total.astype(np.float32),
        v=v_total.astype(np.float32),
        h=h_total.astype(np.float32),
        t=np.array([t], dtype=np.float64),
    )


# -----------------------------
# Build base state and initial perturbation
# -----------------------------
u0, v0, h0, eta0, u0_x, u0_y, v0_x, v0_y = make_rankine_base(X, Y, Umax, R_km * 1000.0, g, H)

# Initial perturbation: Gaussian ring centered at r = R
r = np.sqrt(X**2 + Y**2)
sigma = sigma_km * 1000.0
eta_amp = A_dyn / g  # convert dynamic-height scale to meters

eta_p = eta_amp * np.exp(-((r - R_km * 1000.0)**2) / (2.0 * sigma**2))
u_p = np.zeros_like(eta_p)
v_p = np.zeros_like(eta_p)

# Sponge
sponge = make_sponge(X, Y, Lx, Ly, sponge_width_km * 1000.0, sponge_max)

# Check the requested power radius
Rout = R_OUT_km * 1000.0
if Rout >= 0.5 * min(Lx, Ly):
    raise ValueError(
        f"R_OUT = {R_OUT_km} km does not fit inside a {Lx_km} km square domain centered on the vortex. "
        f"Use R_OUT < {0.5 * min(Lx, Ly) / 1000.0:.1f} km, or enlarge the domain to at least {2*R_OUT_km:.0f} km square."
    )

# -----------------------------
# Time stepping
# -----------------------------
nsteps = int(round(t_end / dt))
output_stride = int(round(3600.0 / dt))  # hourly output
if abs(output_stride * dt - 3600.0) > 1e-12:
    raise ValueError("dt must divide 3600 s exactly for hourly output with this script.")

power_rows = []

# Save t = 0 state
u_total = u0 + u_p
v_total = v0 + v_p
h_total = h0 + eta_p
save_snapshot(0, 0.0, u_total, v_total, h_total)

P0 = compute_power_flux(eta_p, u_p, v_p, Rout, rho, g)
power_rows.append((0.0, P0))

print(f"Grid: {nx} x {ny}, dx = {dx/1000:.1f} km")
print(f"Domain: {Lx_km:.0f} km x {Ly_km:.0f} km")
print(f"Time step: {dt:.1f} s, total steps: {nsteps}")
print(f"Hourly output every {output_stride} steps")
print(f"Power radius: {Rout/1000:.1f} km")
print(f"Initial power flux: {P0:.6e} W")

for n in range(1, nsteps + 1):
    # Strang splitting:
    #   half explicit base-flow coupling
    #   exact gravity-wave step
    #   half explicit base-flow coupling
    eta_p, u_p, v_p = rk4_step(
        eta_p, u_p, v_p, 0.5 * dt,
        lambda ee, uu, vv: rhs_explicit(
            ee, uu, vv,
            u0, v0, eta0,
            u0_x, u0_y, v0_x, v0_y,
            sponge
        )
    )

    eta_p, u_p, v_p = exact_linear_gravity_step(eta_p, u_p, v_p, dt)

    eta_p, u_p, v_p = rk4_step(
        eta_p, u_p, v_p, 0.5 * dt,
        lambda ee, uu, vv: rhs_explicit(
            ee, uu, vv,
            u0, v0, eta0,
            u0_x, u0_y, v0_x, v0_y,
            sponge
        )
    )

    # Optional mild filter to suppress spectral ringing from the Rankine corner and sponge
    # (You can comment this out if you want a strictly unfiltered pseudo-spectral evolution.)
    E = np.fft.fft2(eta_p)
    U = np.fft.fft2(u_p)
    V = np.fft.fft2(v_p)
    kx_cut = (2.0 / 3.0) * np.max(np.abs(kx_1d))
    ky_cut = (2.0 / 3.0) * np.max(np.abs(ky_1d))
    filt = (np.abs(KX) <= kx_cut) & (np.abs(KY) <= ky_cut)
    E *= filt
    U *= filt
    V *= filt
    eta_p = np.fft.ifft2(E).real
    u_p = np.fft.ifft2(U).real
    v_p = np.fft.ifft2(V).real

    t = n * dt

    if n % output_stride == 0 or n == nsteps:
        u_total = u0 + u_p
        v_total = v0 + v_p
        h_total = h0 + eta_p

        save_snapshot(n, t, u_total, v_total, h_total)

        P = compute_power_flux(eta_p, u_p, v_p, Rout, rho, g)
        power_rows.append((t, P))
        print(f"t = {t/3600:.1f} h, power flux = {P:.6e} W")

# Save power time series
power_csv = os.path.join(outdir, "power_flux_timeseries.csv")
with open(power_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["time_s", "power_W"])
    for t, P in power_rows:
        writer.writerow([f"{t:.6f}", f"{P:.12e}"])

print(f"Saved snapshots to: {snapshot_dir}")
print(f"Saved power time series to: {power_csv}")
