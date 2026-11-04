"""Historical 1992 three-layer illustration, not a molecular model.

Reproduces the published modeling assumptions, not experimental data.
The epsilon_2 step is imposed; it does not result from microscopic sampling.
"""
import math
from dataclasses import dataclass
from .molecular import DEBYE, ANGSTROM2, EPSILON_0

@dataclass(frozen=True)
class LegacyPoint:
    area_a2: float
    potential_v: float
    water_v: float
    head_v: float
    tail_v: float
    diffuse_v: float


def legacy_1992_point(area_a2: float) -> LegacyPoint:
    if not 20 <= area_a2 <= 80:
        raise ValueError('Historical parameterization defined for 20 <= A <= 80 A^2')
    a=area_a2
    relative_dielectric_head=6.4 if a<=30 else 78.5
    relative_dielectric_tail=2.8+(a-20.)*(2.18-2.8)/60.
    mu_water_eff=-0.065 * max(0.,min(1.,(30.-a)/10.))
    mu_head_eff=.990 / relative_dielectric_head
    # Published geometry uses a rigid 25 A rod rotating in a cone.
    mu_tail_eff=.330 * math.sqrt(1. - a/(math.pi*25.**2))/relative_dielectric_tail
    factor=DEBYE/(EPSILON_0*ANGSTROM2*a)
    w,h,t=mu_water_eff*factor,mu_head_eff*factor,mu_tail_eff*factor
    psi=-0.0509 * math.asinh(173.0/a)  # published eq. (4), 295 K
    return LegacyPoint(a,w+h+t+psi,w,h,t,psi)
