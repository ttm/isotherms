# Molecular-state theory of fatty-acid Langmuir-film surface potentials

A reproducible **research prototype** of a molecular-level statistical-mechanical explanation for the surface-potential/area isotherm of stearic-acid monolayers.

**Status:** an explicit hypothesis and an illustrative, tested computation. **Not an independently calibrated or experimentally validated model.** There are no raw isotherm datasets in the three source papers uploaded with this task.

The 1992 Oliveira-Taylor-Morgan approach required a step in the effective local headgroup permittivity near 30 Å²/molecule. This project instead represents molecular packing, chain tilt, water orientation and charge dissociation explicitly. The experimental observable is reconstructed from the **excess normal dipole density relative to clean water**, plus a charge-regulated Gouy-Chapman term. No area-dependent dielectric constant is imposed at the headgroup. A separate comparison function reconstructs the historical three-layer approach.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[plot,test]'
pytest -q
python scripts/run_model.py --out results
# Optional: sample correlated 2D molecular domains
python scripts/run_model.py --mc --points 21 --size 32 --sweeps 500 --out results/mc
```

The repository also includes three static example figures under [`figures/`](figures/). The script generates `illustrative_predictions.csv`, `illustrative_isotherm.png`, `molecular_contributions.png`, and `molecular_order.png`. These are **model-only figures**, not fits to published experimental dots.

When an independently measured dataset is available, save it as a CSV with `area_a2,potential_mv,sigma_mv` and run:

```bash
python scripts/fit_data.py my_experimental_isotherm.csv --out results/fit
```

The fit holds most physical parameters fixed, estimates only three free parameters, and reports a local (not globally reliable) covariance. A fit is not evidence of the mechanism without independent measurements of packing or hydration.

## Files

- [`article/molecular_surface_potential.md`](article/molecular_surface_potential.md): research manuscript with derivation, experimental predictions, and limitations.
- [`src/isotherms/molecular.py`](src/isotherms/molecular.py): mean-field equilibrium + 2D checkerboard Metropolis sampler.
- [`src/isotherms/legacy.py`](src/isotherms/legacy.py): reconstruction of the 1992 model for comparison.
- [`scripts/run_model.py`](scripts/run_model.py): curves, CSV, and figures.
- [`scripts/fit_data.py`](scripts/fit_data.py): optional calibration using **real, user-supplied data only**.
- [`tests/test_model.py`](tests/test_model.py): regression, unit, screening and reproducibility checks.

## Scope and epistemic status

The lattice sites represent individual fatty-acid molecules and one effective hydration-water orientational degree of freedom each. The model is *molecular-state coarse-graining*, not force-field atomistic molecular dynamics, and its water parameters are illustrative, not measured. Bulk water permittivity **is still used**, appropriately, in the diffuse Gouy-Chapman approximation. An unqualified claim that this is a full microscopic electrostatic solution would be incorrect.

The 2005 de Moura, Oliveira Jr., Trsic molecular-dynamics abstract already addressed molecular contributions and reported a roughly six-fold overestimate of measured potential. Novelty here is a mechanistic, area-dependent, charge-regulated molecular ensemble and experimentally discriminating predictions, **not** the first molecular investigation of Langmuir-film surface potential.

## Main references

1. Oliveira, Taylor & Morgan (1992), *Thin Solid Films* **210/211**, 76-78. https://doi.org/10.1016/0040-6090(92)90172-8
2. Taylor, Oliveira & Morgan (1989), *Chemical Physics Letters* **161**, 147-150.
3. Taylor, Oliveira & Morgan (1990), *Journal of Colloid and Interface Science* **139**, 508-518.
4. Demchak & Fort (1974), *J. Colloid Interface Sci.* **46**, 191-202. https://doi.org/10.1016/0021-9797(74)90002-2
5. Taylor & Bayes (1994), *Physical Review E* **49**, 1439. https://doi.org/10.1103/PhysRevE.49.1439
6. de Moura, Oliveira Jr. & Trsic (2005), abstract in XIII SBQT proceedings. https://marcia.iqm.unicamp.br/SBQT13.pdf
7. Molecular-dynamics study of hydrated DPPC monolayers with planar charge-density integration: https://pmc.ncbi.nlm.nih.gov/articles/PMC3488352/

Re-use and redistribution of the three user-supplied copyrighted articles is intentionally avoided. Only original exposition and citations are committed.
