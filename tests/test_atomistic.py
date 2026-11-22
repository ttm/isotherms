import numpy as np
import pytest

from isotherms.atomistic import (
    DEBYE, ELEMENTARY_CHARGE, EPSILON_0, NANOMETER,
    planar_poisson, surface_step, excess_surface_potential, molecular_area_a2,
)


def symmetrical_slab(q=0.012, separation=0.3, area_a2=25.0, n=64, bins=2048):
    """Two symmetric outward dipole sheets: neutral 3D-periodic slab benchmark."""
    Lxy = np.sqrt(n * area_a2 / 100.0)  # nm
    box = (Lxy, Lxy, 14.0)
    # Upper outward dipoles at z=9, lower mirrored outward dipoles at z=3.
    # Lower: positive at 3-separation/2, negative at 3+separation/2.
    # Upper: negative at 9-separation/2, positive at 9+separation/2.
    z = np.concatenate([
        np.full(n, 3-separation/2), np.full(n, 3+separation/2),
        np.full(n, 9-separation/2), np.full(n, 9+separation/2),
    ])
    charges = np.concatenate([
        np.full(n, q), np.full(n, -q), np.full(n, -q), np.full(n, q),
    ])
    return planar_poisson(z, charges, box_nm=box, bins=bins), box


def test_zero_reference():
    p = planar_poisson([1.0, 1.0], [0.5, -0.5], box_nm=(5, 5, 14))
    assert np.max(np.abs(p.potential_v)) < 1e-12
    assert np.max(np.abs(p.charge_density_c_m3)) < 1e-3


def test_charge_neutrality_enforced():
    with pytest.raises(ValueError, match="Nonneutral"):
        planar_poisson([2, 3], [1, -0.7], box_nm=(5,5,14))


def test_helmholtz_dipole_sheet_analytic():
    profile, box = symmetrical_slab()
    voltage = surface_step(profile, bulk_water_nm=(5,7), adjacent_air_nm=(10,11))
    area_m2 = box[0]*box[1]*NANOMETER**2
    # Exact infinite sheet plateau step: n*q*e*d/(eps0*S)
    exact = (64*0.012*ELEMENTARY_CHARGE * 0.3*NANOMETER)/(EPSILON_0*area_m2)
    assert voltage == pytest.approx(exact, rel=0.015)
    clean = planar_poisson([1, 1], [0.0, 0.0], box_nm=box, bins=2048)
    assert excess_surface_potential(profile, clean, bulk_water_nm=(5,7), adjacent_air_nm=(10,11)) == pytest.approx(exact, rel=0.015)


def test_area_inverse_scaling():
    v=[]
    for a in (20, 40, 60):
        p, box = symmetrical_slab(area_a2=a)
        assert molecular_area_a2(box, 64) == pytest.approx(a)
        v.append(surface_step(p, bulk_water_nm=(5,7), adjacent_air_nm=(10,11)))
    assert np.array(v) / v[0] == pytest.approx([1, .5, 1/3], rel=0.02)


def test_potential_invariant_to_charge_order_and_periodic_wrap():
    p, box = symmetrical_slab()
    z = np.array([9.15, 8.85, 3.15, 2.85])
    q = np.array([0.5, -0.5, -0.5, 0.5])
    pa = planar_poisson(z,q,box_nm=box)
    pb = planar_poisson(z[[2,0,3,1]]+14, q[[2,0,3,1]], box_nm=box)
    assert np.allclose(pa.potential_v, pb.potential_v, atol=1e-12)


def test_water_background_subtraction():
    p, box = symmetrical_slab()
    # Add a charge-neutral background with the same box and atom model in both,
    # but here the clean array is used to illustrate that only the difference matters.
    background = planar_poisson([3.7,3.95,8.05,8.3],[.02,-.02,.02,-.02],box_nm=box)
    clean = background
    film_bg = planar_poisson(
        np.r_[np.repeat([2.85, 3.15, 8.85, 9.15],64), [3.7,3.95,8.05,8.3]],
        np.r_[np.repeat([.012,-.012,-.012,.012],64), [.02,-.02,.02,-.02]],
        box_nm=box
    )
    v=excess_surface_potential(film_bg,clean,bulk_water_nm=(5,7),adjacent_air_nm=(10,11))
    expected=surface_step(p,bulk_water_nm=(5,7),adjacent_air_nm=(10,11))
    assert v == pytest.approx(expected,rel=1e-8)


def test_excess_rejects_mismatched_cell():
    p,_=symmetrical_slab()
    q=planar_poisson([2,2],[1,-1],box_nm=(2,2,12))
    with pytest.raises(ValueError,match="same z grid"):
        excess_surface_potential(p,q,bulk_water_nm=(5,7),adjacent_air_nm=(10,11))


def test_framewise_charge_conservation_and_variable_charge_frames():
    z = np.array([[3.0,3.2,9.0,9.2], [3.0,3.2,9.0,9.2]])
    q = np.array([[.01,-.01,-.01,.01], [.02,-.02,-.02,.02]])
    v_avg = planar_poisson(z,q,box_nm=(5,5,14))
    v1 = planar_poisson(z[0],q[0],box_nm=(5,5,14))
    v2 = planar_poisson(z[1],q[1],box_nm=(5,5,14))
    assert np.allclose(v_avg.potential_v, (v1.potential_v+v2.potential_v)/2)
    q[1,0] += .1
    with pytest.raises(ValueError,match="Nonneutral"):
        planar_poisson(z,q,box_nm=(5,5,14))
