import numpy as np
import pytest
from dataclasses import replace
from isotherms.molecular import (Parameters, mean_field_point, metropolis_point,
                                charge_regulation, diffuse_potential, _tilt_cos)


def test_monotonic_packing_and_physical_bounds():
    areas = [20, 25, 30, 35, 40, 60, 80]
    points = [mean_field_point(a) for a in areas]
    fractions = [v.condensed_fraction for v in points]
    assert all(x >= y for x,y in zip(fractions,fractions[1:]))
    assert fractions[0] > 0.98
    assert fractions[-1] < 0.02
    for v in points:
        assert 0 <= v.dissociation <= 1
        assert -1 <= v.water_order <= 1
        assert 0 <= v.tail_cosine <= 1
        assert v.double_layer_v < 0
        assert v.potential_v == pytest.approx(v.head_v+v.tail_v+v.water_v+v.double_layer_v)


def test_charge_regulation_and_screening():
    p=Parameters()
    alpha,psi=charge_regulation(20,.9,-.4,p)
    assert alpha == pytest.approx(0.0021876473508838346,rel=0.001)
    assert psi == pytest.approx(-.14979351923234133,rel=0.001)
    p_salt=replace(p,ionic_strength_m=.01)
    alpha2,psi2=charge_regulation(20,.9,-.4,p_salt)
    assert alpha2 > alpha
    assert abs(psi2) < abs(psi)
    assert diffuse_potential(0,20,p) == 0


def test_water_coupling_affects_potential():
    p=Parameters()
    no_water=replace(p,hydration_water_d=0.)
    assert mean_field_point(20,no_water).potential_v > mean_field_point(20,p).potential_v
    assert _tilt_cos(.5,20.) == pytest.approx(.5,abs=.002)


def test_qualitative_isotherm_not_claimed_fit():
    high=mean_field_point(80)
    mid=mean_field_point(40)
    low=mean_field_point(20)
    # Broad, physically grounded diagnostics; NOT comparison to digitized experiments.
    assert abs(high.potential_v) < .10
    assert abs(mid.potential_v) < .15
    assert low.potential_v > mid.potential_v+.10


def test_mc_is_finite_reproducible():
    a,sa,wa=metropolis_point(30,size=16,burn_sweeps=40,sample_sweeps=40,seed=12)
    b,sb,wb=metropolis_point(30,size=16,burn_sweeps=40,sample_sweeps=40,seed=12)
    assert a == b
    assert np.array_equal(sa,sb)
    assert np.array_equal(wa,wb)
    assert set(np.unique(sa)).issubset({0,1})
    assert set(np.unique(wa)).issubset({-1,1})


def test_historical_model_expanded_and_condensed():
    from isotherms.legacy import legacy_1992_point
    assert legacy_1992_point(20).potential_v > .15
    assert abs(legacy_1992_point(80).potential_v) < .15
    with pytest.raises(ValueError):legacy_1992_point(90)


def test_fit_data_validation(tmp_path):
    from scripts.fit_data import load
    data=tmp_path/'bad.csv'
    data.write_text('area_a2,potential_mv,sigma_mv\n20,100,10\n')
    with pytest.raises(ValueError):load(data)
