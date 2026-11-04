"""Coupled molecular-state / interfacial-water model for fatty-acid Langmuir films.

Each site represents one fatty acid and one *excess* first-shell water orientation.
Two molecular packing states (expanded/condensed) interact laterally; the
water direction is an explicit two-state molecular orientation. The tail
orientation is thermally distributed. Acid ionization obeys charge-regulated
Gouy-Chapman theory. No area-dependent effective permittivity is imposed.

This is a coarse-grained molecular model, not an atomistic MD simulation.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar

EPSILON_0 = 8.8541878128e-12  # F / m
ELEMENTARY_CHARGE = 1.602176634e-19  # C
BOLTZMANN = 1.380649e-23  # J / K
AVOGADRO = 6.02214076e23  # mol^-1
DEBYE = 3.33564e-30  # C m
ANGSTROM2 = 1e-20  # m^2

@dataclass(frozen=True)
class Parameters:
    """All couplings are in k_B T; all dipoles are in Debye.

    Illustrative values, *not fitted or independently determined*.
    Mean-field is single-phase at j_chain < 1/4 for square-lattice geometry.
    """
    temperature_k: float = 295.0
    ionic_strength_m: float = 2.5e-6  # effective monovalent electrolyte ionic strength
    water_permittivity: float = 78.5  # bulk permittivity for diffuse PB *only*
    bulk_ph: float = 5.6
    pka0: float = 5.7
    pka_order_shift: float = 0.0  # unverified coupling kept optional
    pka_water_shift: float = 0.0
    critical_area_a2: float = 30.0
    packing_slope_per_a2: float = 0.45
    j_chain: float = 0.20  # Ising sigma*sigma NN, dimensionless
    j_water: float = 0.0  # enable in MC only; MF uses 0 by assumption
    water_bias_at_80: float = -0.37
    water_overlap_bias: float = -0.12
    water_condensed_shift: float = 0.09
    head_protonated_d: float = 0.99
    head_ionized_d: float = 0.75
    tail_group_d: float = 0.33
    hydration_water_d: float = 2.5  # effective excess dipole / first shell water
    expanded_tilt_cos_target: float = 0.48
    condensed_tilt_cos_target: float = 1.02
    expanded_tilt_stiffness: float = 8.0
    condensed_tilt_stiffness: float = 150.0

    def validate(self) -> None:
        if self.temperature_k <= 0 or self.ionic_strength_m <= 0 or self.water_permittivity <= 0:
            raise ValueError('T, ionic strength, and bulk permittivity must be positive')
        if self.critical_area_a2 <= 0 or self.packing_slope_per_a2 <= 0:
            raise ValueError('critical area and packing slope must be positive')
        if self.expanded_tilt_stiffness <= 0 or self.condensed_tilt_stiffness <= 0:
            raise ValueError('tilt stiffnesses must be positive')
        if self.j_water != 0:
            raise ValueError('mean-field curve requires j_water=0; MC permits nonzero j_water')


def water_fields(area_a2: float, p: Parameters) -> tuple[float, float]:
    """The second field is for condensed molecules; bias is not dielectric epsilon."""
    if area_a2 <= 0:
        raise ValueError('Area per molecule must be positive')
    be = p.water_bias_at_80 + p.water_overlap_bias * (40.0 / area_a2 - 0.5)
    return be, be + p.water_condensed_shift


def packing_bias(area_a2: float, p: Parameters) -> float:
    return p.packing_slope_per_a2 * (area_a2 - p.critical_area_a2)


def _xlogx(x: float) -> float:
    return 0.0 if x <= 0.0 else x * math.log(x)


def _mean_field_order(area_a2: float, p: Parameters) -> float:
    be, bc = water_fields(area_a2, p)
    # Exact single-site trace over water orientations (j_water == 0).
    effective_bias = packing_bias(area_a2, p) - (np.logaddexp(bc, -bc) - np.logaddexp(be, -be))
    def f(m: float) -> float:
        pc = (1.0 + m) / 2.0
        pe = 1.0 - pc
        return (pc * effective_bias - 2.0 * p.j_chain * m * m
                + _xlogx(pc) + _xlogx(pe))
    # Minimize all local basins, not a single fixed-point iteration, so the
    # globally stable branch is chosen even for first-order transitions.
    grid = np.linspace(-1.0, 1.0, 1001)
    vals = np.array([f(float(x)) for x in grid])
    ids = np.r_[0, np.where((vals[1:-1] <= vals[:-2]) & (vals[1:-1] <= vals[2:]))[0]+1, len(grid)-1]
    candidates = [(float(vals[i]),float(grid[i])) for i in ids]
    for i in ids:
        if 0 < i < len(grid)-1:
            res = minimize_scalar(f, method='bounded', bounds=(grid[i-1],grid[i+1]),
                                  options={'xatol':1e-12})
            candidates.append((float(res.fun),float(res.x)))
    return min(candidates)[1]


@lru_cache(maxsize=64)
def _tilt_cos(target: float, stiffness: float) -> float:
    # Continuous molecular tilt degrees of freedom u = cos(theta), 0 <= u <= 1.
    nodes, weights = np.polynomial.legendre.leggauss(96)
    u = 0.5 * (nodes + 1.)
    logw = -0.5*stiffness*(u-target)**2
    density = np.exp(logw - np.max(logw)) * weights
    return float(np.sum(u * density) / np.sum(density))


def diffuse_potential(alpha: float, area_a2: float, p: Parameters) -> float:
    """Grahame relation, flat 1:1 Gouy-Chapman layer; negative for A- headgroups."""
    if not (0.0 <= alpha <= 1.0) or area_a2 <= 0:
        raise ValueError('alpha must be between 0 and 1, area > 0')
    nb = p.ionic_strength_m * 1000.0 * AVOGADRO
    sigma = -ELEMENTARY_CHARGE * alpha / (area_a2*ANGSTROM2)
    charge_scale = math.sqrt(8.0 * EPSILON_0 * p.water_permittivity * BOLTZMANN * p.temperature_k * nb)
    return (2.0 * BOLTZMANN * p.temperature_k / ELEMENTARY_CHARGE) * math.asinh(sigma / charge_scale)


def charge_regulation(area_a2: float, phi: float, water_order: float, p: Parameters) -> tuple[float, float]:
    """Solve Henderson-Hasselbalch with local proton activity at psi0.

    logit(alpha) = ln(10)*(pH - pKa) + beta*e*psi0.
    """
    pka = p.pka0 + p.pka_order_shift*phi + p.pka_water_shift*water_order
    def residual(alpha: float) -> float:
        psi = diffuse_potential(alpha, area_a2, p)
        return math.log(alpha/(1.-alpha)) - math.log(10.)*(p.bulk_ph-pka) - (ELEMENTARY_CHARGE*psi/(BOLTZMANN*p.temperature_k))
    alpha = brentq(residual, 1e-13, 1.-1e-12, xtol=1e-13)
    return alpha, diffuse_potential(alpha, area_a2, p)


@dataclass(frozen=True)
class Point:
    area_a2: float
    condensed_fraction: float
    water_order: float
    tail_cosine: float
    dissociation: float
    double_layer_v: float
    head_v: float
    tail_v: float
    water_v: float
    potential_v: float


def _build_point(area_a2: float, phi: float, water_order: float, p: Parameters) -> Point:
    if not (0 <= phi <= 1):
        raise ValueError('condensed fraction out of range')
    tail_cos = ((1-phi)*_tilt_cos(p.expanded_tilt_cos_target,p.expanded_tilt_stiffness)
                + phi*_tilt_cos(p.condensed_tilt_cos_target,p.condensed_tilt_stiffness))
    alpha, psi = charge_regulation(area_a2, phi, water_order, p)
    factor = DEBYE / (EPSILON_0 * ANGSTROM2 * area_a2)
    head = ((1-alpha)*p.head_protonated_d + alpha*p.head_ionized_d)*factor
    tail = p.tail_group_d*tail_cos*factor
    water = p.hydration_water_d*water_order*factor
    return Point(area_a2,phi,water_order,tail_cos,alpha,psi,head,tail,water,head+tail+water+psi)


def mean_field_point(area_a2: float, p: Parameters=Parameters()) -> Point:
    """Global equilibrium mean-field solution with exact local water trace."""
    p.validate()
    m = _mean_field_order(area_a2,p)
    phi = (1.+m)/2.
    be,bc = water_fields(area_a2,p)
    water = (1.-phi)*math.tanh(be) + phi*math.tanh(bc)
    return _build_point(area_a2,phi,water,p)


def mean_field_curve(areas_a2: np.ndarray, p: Parameters=Parameters()) -> list[Point]:
    return [mean_field_point(float(a),p) for a in areas_a2]


def _neighbours(arr: np.ndarray) -> np.ndarray:
    return np.roll(arr,1,axis=0)+np.roll(arr,-1,axis=0)+np.roll(arr,1,axis=1)+np.roll(arr,-1,axis=1)


def metropolis_point(area_a2: float, p: Parameters=Parameters(), size: int=32,
                     burn_sweeps: int=300, sample_sweeps: int=300, seed: int=2026) -> tuple[Point,np.ndarray,np.ndarray]:
    """Two-dimensional checkerboard Metropolis; samples correlated molecular domains.

    Returns (mean macroscopic point, final packing map, final water-orientation map).
    This is a finite lattice; uncertainty needs independent seeds, not just sweeps.
    """
    # j_water is allowed in MC, unlike the exactly traced mean-field estimator.
    if area_a2 <= 0 or size < 4 or size%2 or burn_sweeps < 0 or sample_sweeps < 1:
        raise ValueError('size must be even >= 4, area > 0, samples >= 1')
    rng = np.random.default_rng(seed)
    s = rng.integers(0,2,(size,size),dtype=np.int8)
    w = rng.choice(np.array([-1,1],dtype=np.int8),size=(size,size))
    x,y=np.indices((size,size))
    masks=((x+y)%2==0,(x+y)%2==1)
    be,bc=water_fields(area_a2,p)
    d=packing_bias(area_a2,p)
    sums=np.zeros(2,dtype=float)
    for sweep in range(burn_sweeps+sample_sweeps):
        for mask in masks:
            sigma=2*s-1
            ds=(1-2*s)*(d-2*p.j_chain*_neighbours(sigma)-p.water_condensed_shift*w)
            flip=mask & ((ds<=0)|(rng.random((size,size))<np.exp(-np.maximum(ds,0))))
            s[flip]=1-s[flip]
            dw=2*w*(be+p.water_condensed_shift*s+p.j_water*_neighbours(w))
            flip=mask & ((dw<=0)|(rng.random((size,size))<np.exp(-np.maximum(dw,0))))
            w[flip]=-w[flip]
        if sweep>=burn_sweeps:
            sums += [np.mean(s),np.mean(w)]
    phi,water=sums/sample_sweeps
    return _build_point(area_a2,float(phi),float(water),p),s,w
