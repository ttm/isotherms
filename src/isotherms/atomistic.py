"""Atom-resolved, parameter-free (at analysis stage) electrostatic forward operator.

Consumes equilibrated and independently parameterized atomistic coordinates and charges.
It does NOT generate molecular dynamics, validate a force field, or fit an isotherm.

Periodic planar Poisson boundary condition: v(k) = rho(k)/(eps0*k*k), k != 0.
Use a symmetric two-sided, neutral slab, with adequately wide bulk-water and
vacuum plateaus; asymmetric slab electrostatics require an explicit correction
consistent with the simulation engine. No dielectric fitting is performed.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

EPSILON_0 = 8.8541878128e-12  # F/m
ELEMENTARY_CHARGE = 1.602176634e-19  # C
NANOMETER = 1e-9
ANGSTROM_SQUARED = 1e-20
DEBYE = 3.33564e-30  # C m


@dataclass(frozen=True)
class ChargeProfile:
    z_nm: np.ndarray
    charge_density_c_m3: np.ndarray
    potential_v: np.ndarray
    electric_field_v_m: np.ndarray
    net_charge_e: float


def planar_poisson(
    z_nm: np.ndarray,
    charges_e: np.ndarray,
    *,
    box_nm: tuple[float, float, float],
    bins: int = 1024,
    neutral_tolerance_e: float = 1e-6,
) -> ChargeProfile:
    """CIC-deposit point charges and solve the periodic, one-dimensional Poisson equation.

    ``z_nm`` may be N atoms or (frames, N atoms); ``charges_e`` is
    (atoms,) for fixed charges or (frames, atoms) for variable charges.
    The charge density and electrostatic potential are frame-averaged. Charge
    is neutral in EACH frame. A constant potential is immaterial; differences
    between well-separated plateaus have a physical meaning for matching BC.
    """
    z = np.asarray(z_nm, dtype=np.float64)
    q = np.asarray(charges_e, dtype=np.float64)
    if z.ndim == 1:
        z = z[None, :]
    if q.ndim == 1:
        q = np.broadcast_to(q, z.shape) if z.ndim == 2 and q.size == z.shape[1] else q
    if z.ndim != 2 or q.shape != z.shape or z.shape[0] == 0:
        raise ValueError("z_nm and charges_e must have compatible (frames, atoms) shapes")
    Lx, Ly, Lz = map(float, box_nm)
    if not all(np.isfinite(x) and x > 0 for x in (Lx, Ly, Lz)):
        raise ValueError("box_nm must have positive finite dimensions")
    if bins < 16 or int(bins) != bins:
        raise ValueError("bins must be an integer >= 16")
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(q)):
        raise ValueError("charges and positions must be finite")
    frame_nets = np.sum(q, axis=1)
    if np.any(np.abs(frame_nets) > neutral_tolerance_e):
        raise ValueError(f"Nonneutral cell (max |Q|={np.max(np.abs(frame_nets)):.6g} e); add explicit counterions")
    net = float(np.mean(frame_nets))

    dz_nm = Lz / bins
    raw = np.mod(z, Lz) / dz_nm
    left = np.floor(raw).astype(np.int64) % bins
    right = (left + 1) % bins
    fraction = raw - np.floor(raw)
    mass = np.zeros(bins, dtype=np.float64)
    for f in range(z.shape[0]):
        np.add.at(mass, left[f], q[f] * (1.0 - fraction[f]))
        np.add.at(mass, right[f], q[f] * fraction[f])
    mass /= z.shape[0]
    voxel_m3 = (Lx * Ly * dz_nm) * NANOMETER**3
    rho = mass * ELEMENTARY_CHARGE / voxel_m3
    # Remove round-off only, not physical nonneutrality.
    rho -= rho.mean()
    dz_m = dz_nm * NANOMETER
    k = 2 * np.pi * np.fft.rfftfreq(bins, d=dz_m)
    rho_hat = np.fft.rfft(rho)
    phi_hat = np.zeros_like(rho_hat)
    phi_hat[1:] = rho_hat[1:] / (EPSILON_0 * k[1:]**2)
    phi = np.fft.irfft(phi_hat, n=bins)
    e_hat = np.zeros_like(rho_hat)
    e_hat[1:] = -1j * k[1:] * phi_hat[1:]
    e_field = np.fft.irfft(e_hat, n=bins)
    z_axis_nm = np.arange(bins, dtype=float) * dz_nm
    return ChargeProfile(z_axis_nm, rho, phi, e_field, net)


def plateau_mean(profile: ChargeProfile, window_nm: tuple[float, float]) -> float:
    a, b = map(float, window_nm)
    z = profile.z_nm
    if not (0 <= a < b <= z.size * (z[1] - z[0])):
        raise ValueError("plateau window outside box or reversed")
    mask = (z >= a) & (z < b)
    if mask.sum() < 3:
        raise ValueError("plateau window too narrow for chosen bins")
    return float(np.mean(profile.potential_v[mask]))


def surface_step(
    profile: ChargeProfile,
    *,
    bulk_water_nm: tuple[float, float],
    adjacent_air_nm: tuple[float, float],
) -> float:
    """Potential in adjacent air minus potential in central bulk water (volts)."""
    return plateau_mean(profile, adjacent_air_nm) - plateau_mean(profile, bulk_water_nm)


def excess_surface_potential(
    film: ChargeProfile,
    clean_water: ChargeProfile,
    *,
    bulk_water_nm: tuple[float, float],
    adjacent_air_nm: tuple[float, float],
) -> float:
    """Monolayer minus clean-water step, same electrostatic convention and box.

    The caller must ensure identical binning, dimensions and consistent atomistic
    force-field / charge model. This never infers a local dielectric constant.
    """
    if film.z_nm.shape != clean_water.z_nm.shape or not np.allclose(
        film.z_nm, clean_water.z_nm, rtol=0, atol=1e-12
    ):
        raise ValueError("film and clean profiles must use the same z grid")
    return surface_step(film, bulk_water_nm=bulk_water_nm, adjacent_air_nm=adjacent_air_nm) - surface_step(
        clean_water, bulk_water_nm=bulk_water_nm, adjacent_air_nm=adjacent_air_nm
    )


def molecular_area_a2(box_nm: tuple[float, float, float], n_per_interface: int) -> float:
    if not isinstance(n_per_interface, int) or n_per_interface <= 0:
        raise ValueError("n_per_interface must be a positive integer")
    return box_nm[0] * box_nm[1] * 100.0 / n_per_interface
