"""Check provenance, numerical integrity and determinism of the comparison."""
from __future__ import annotations
from pathlib import Path
import importlib.util
import numpy as np
import pytest

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/compare_historical_isotherm.py'
spec=importlib.util.spec_from_file_location('historical_comparison',SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def _tables():
    return mod.read_table(mod.EXPERIMENT),mod.read_table(mod.PILOT)


def test_trace_is_explicitly_approximate_and_unfitted():
    exp,pilot=_tables()
    r=mod.compare(exp,pilot)
    assert r['valid_model_prediction'] is False
    assert r['conditions_match'] is False
    assert 'NOT raw' in r['experimental_provenance']
    assert len(r['comparison'])==6
    assert list(sorted(exp['area_a2']))==list(exp['area_a2'])
    assert exp['potential_mv'][0] > 200
    assert abs(exp['potential_mv'][-1]) < 20
    assert max(pt['raw_difference_mv'] for pt in r['comparison'])>900


def test_interpolation_does_not_fit_or_rezero_model():
    exp,pilot=_tables()
    comparison=mod.compare(exp,pilot)['comparison']
    for i,row in enumerate(comparison):
        assert row['pilot_mv'] == pytest.approx(1000*pilot['voltage_average_v'][i])
        assert row['raw_difference_mv'] == pytest.approx(row['pilot_mv']-np.interp(
            row['area_a2'],exp['area_a2'],exp['potential_mv']))


def test_reject_bad_trace():
    exp,pilot=_tables()
    broken=dict(exp)
    broken['area_a2']=exp['area_a2'].copy()
    broken['area_a2'][2]=broken['area_a2'][1]
    with pytest.raises(ValueError):mod.compare(broken,pilot)


def test_figure_and_diagnostics_created(tmp_path):
    exp,pilot=_tables()
    png=tmp_path/'plot.png'
    r=mod.plot(exp,pilot,png)
    assert png.exists() and png.stat().st_size>10000
    assert not r['valid_model_prediction']