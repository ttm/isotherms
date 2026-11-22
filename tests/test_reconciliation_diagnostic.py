"""Numerical integrity tests, not validation of a molecular mechanism."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "reconciliation_diagnostic.py"
spec = importlib.util.spec_from_file_location("reconciliation_diagnostic", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_calibration_diagnostic_from_repository_data():
    historical=module._table(module.EXPERIMENT)
    md=module._table(module.MD)
    report, arrays=module.diagnose(historical,md)
    assert report['validated_molecular_prediction'] is False
    assert report['conditions_match'] is False
    assert report['n_md_areas']==6
    assert report['raw_rmse_mv'] > report['fit_affine_rmse_mv']
    assert report['fit_affine_rmse_mv'] == pytest.approx(44.934, rel=0.002)
    assert report['fit_affine_gain'] == pytest.approx(0.2590451, rel=1e-5)
    assert report['fit_affine_offset_mv'] == pytest.approx(-132.555, rel=2e-5)
    assert report['fit_affine_loocv_rmse_mv'] > report['fit_affine_rmse_mv']
    assert any((a,b)==(30.,32.) for a,b,_,_ in report['opposite_slope_intervals_a2'])
    assert np.ptp(arrays['md_raw_mv']) > 1000
    assert not np.array_equal(arrays['md_raw_mv'], arrays['affine_ols_mv'])


def test_gaussian_area_smoothing_is_descriptive():
    area=np.array([24.,28.,30.,32.,36.,45.])
    series=np.array([10.,8.,6.,9.,4.,2.])
    smoothed=module.gaussian_area_average(area,series,2.)
    assert smoothed.shape==series.shape
    assert np.all(np.isfinite(smoothed))
    assert smoothed.max() <= series.max()
    assert smoothed.min() >= series.min()
    with pytest.raises(ValueError):
        module.gaussian_area_average(area,series,-1.)


def test_invalid_trace_rejected():
    e=dict(area_a2=np.array([20.,30.,40.,50.,60.]),potential_mv=np.arange(5.),
           reading_uncertainty_mv=np.ones(5))
    m=dict(area_a2=np.array([24.,30.,36.,45.]),voltage_average_v=np.ones(4))
    bad=dict(e,reading_uncertainty_mv=-np.ones(5))
    with pytest.raises(ValueError):module.diagnose(bad,m)
    bad=dict(e,area_a2=e['area_a2'][::-1])
    with pytest.raises(ValueError):module.diagnose(bad,m)