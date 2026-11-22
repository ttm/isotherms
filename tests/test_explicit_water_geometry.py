"""Geometry invariants for explicit-water pilot; do not require OpenMM/forcefield."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip('rdkit')
SCRIPT=Path(__file__).resolve().parents[1]/'scripts/explicit_water_pilot.py'
spec=importlib.util.spec_from_file_location('explicit_water_pilot',SCRIPT)
pilot=importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)

@pytest.mark.parametrize('area',[24,28,30,32,36,45])
def test_stearic_geometry_and_density(area):
    g=pilot.make_coordinates(area,side=4,seed=41)
    assert g['n_acid_per_leaflet']==16
    assert g['n_acid_atoms']==32*56
    assert np.isclose(g['box_nm'][0]*g['box_nm'][1]*100/16,area)
    assert g['n_water']>100
    assert len(g['positions_nm'])==g['n_acid_atoms']+3*g['n_water']
    assert len(g['elements'])==len(g['positions_nm'])
    assert len(g['residues'])==32+g['n_water']
    assert np.all(np.isfinite(g['positions_nm']))
    assert np.ptp(g['positions_nm'][:56,2])>2.0  # extended C18 start
    assert np.all(g['positions_nm']>=-0.5)
    assert np.all(g['positions_nm'][:,2]<g['box_nm'][2])


def test_reproducibility_and_fixed_protonation():
    a=pilot.make_coordinates(30,seed=41)
    b=pilot.make_coordinates(30,seed=41)
    np.testing.assert_array_equal(a['positions_nm'],b['positions_nm'])
    assert len([r for r in a['residues'] if r[0]=='SA'])==32
    assert pilot.Chem.rdMolDescriptors.CalcMolFormula(a['mol'])=='C18H36O2'


def test_invalid_geometry():
    for kwargs in [dict(area_a2=12),dict(area_a2=30,side=1),dict(area_a2=30,box_z_nm=7)]:
        with pytest.raises(ValueError):pilot.make_coordinates(**kwargs)