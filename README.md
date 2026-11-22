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

The script generates `illustrative_predictions.csv`, `illustrative_isotherm.png`, `molecular_contributions.png`, and `molecular_order.png`. These are **model-only figures**, not fits to published experimental dots.

When an independently measured dataset is available, save it as a CSV with `area_a2,potential_mv,sigma_mv` and run:

```bash
python scripts/fit_data.py my_experimental_isotherm.csv --out results/fit
```

The fit holds most physical parameters fixed, estimates only three free parameters, and reports a local (not globally reliable) covariance. A fit is not evidence of the mechanism without independent measurements of packing or hydration.

## Files

- [`article/molecular_surface_potential.md`](article/molecular_surface_potential.md): research manuscript with derivation, experimental predictions, and limitations.
- [`article/molecular_surface_potential.pdf`](article/molecular_surface_potential.pdf): typeset article PDF, including equations, tables and the three computed figures.
- [`scripts/build_article_pdf.sh`](scripts/build_article_pdf.sh): rebuild the PDF from the Markdown source (requires Pandoc, XeLaTeX and Noto Serif).
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


## Atom-resolved forward theory (22 November 2026)

The second research manuscript pursues a stronger, **prediction-first** treatment: derive a surface-potential isotherm from independently specified atomistic charge distributions of stearic acid, explicit water, and ions, rather than fitting an area-dependent headgroup dielectric constant or prescribing an onset near 30 Å².

- [New article, Markdown](article/atomistic_surface_potential.md) and [PDF](article/atomistic_surface_potential.pdf)
- [Planar charge-density/Poisson solver](src/isotherms/atomistic.py), [trajectory analysis CLI](scripts/atomistic_forward.py), [eight regression tests](tests/test_atomistic.py)
- [Analytic sheet benchmark](figures/atomistic_sheet_benchmark.png) and [CSV](examples/atomistic_sheet_benchmark.csv)
- PDF rebuild: `bash scripts/build_atomistic_pdf.sh` (requires Pandoc + pdfLaTeX).

```bash
PYTHONPATH=src pytest -q tests/test_atomistic.py
PYTHONPATH=src python scripts/atomistic_forward.py --demo --out results/atomistic
# For real equilibrated, charge-neutral atomistic trajectories:
PYTHONPATH=src python scripts/atomistic_forward.py --film film_A30.npz --clean clean_A30.npz --water 5 7 --air 10 11
```

**Validation scope:** electrostatic synthetic-sheet and neutrality tests pass. A later nonconverged stearic-acid/water MD pilot ran successfully, but cannot establish an equilibrium isotherm, predictive onset, or agreement with measurements. See the quantitative quality audit below. The original molecular-state article and model remain intact.

## Explicit-water MD pilot and scientific QA (22 November 2026)

The [six-area OpenMM pilot](scripts/explicit_water_pilot.py) executes neutral stearic-acid monolayers with explicit TIP3P water at 24, 28, 30, 32, 36, and 45 Å²/molecule. It ran successfully in [GitHub Actions](https://github.com/ttm/isotherms/actions/runs/45057869556), with [`summary.csv`](results/explicit_water_pilot/summary.csv), [`protocol.json`](results/explicit_water_pilot/protocol.json), and a [pilot curve](results/explicit_water_pilot/pilot_isotherm.png). Full trajectory NPZ files and charge profiles are retained as an Actions artifact for 30 days.

**This is a nonconverged execution check, not an atomistically predicted surface-potential isotherm.** The measured numerical results have substantial frame fluctuations, disagreements between nominally equivalent interfaces, and non-flat nominal bulk-water reference windows. A quantitative [scientific quality audit](article/explicit_water_pilot_quality_audit.md) reports these failures with numbers extracted from the archived trajectories. See [molecular modeling alternatives](article/molecular_modeling_routes.md) for cDFT, quantum-derived lattice statistics, polarizable MD and other routes.

Workflow summary: [regular tests](https://github.com/ttm/isotherms/actions/runs/45057869567) passed **15 tests, with one optional RDKit geometry module skipped** in that environment; the OpenMM workflow independently passed all **8 geometry tests** before running six simulations. Green CI means the software executed, **not** that the predicted voltage is scientifically valid. No experimental fit or independent physical validation is claimed.

To repeat locally, use the [micromamba dependency setup](.github/workflows/explicit-water-pilot.yml) and run `PYTHONPATH=src python scripts/explicit_water_pilot.py --areas 24 28 30 32 36 45 --side 4 --equil-ps 5 --prod-ps 10`. This is computationally expensive relative to the unit tests and currently uses unpinned binary package versions; archive the resolved environment before production research.
## Historical surface-potential benchmark (22 November 2026)

The [atom-resolved article](article/atomistic_surface_potential.pdf) now includes
[a direct, exploratory graphical comparison](figures/experimental_vs_md_pilot.png)
against Figure 1 of Oliveira, Taylor & Morgan (1992). The reference
[data CSV](data/oliveira_1992_fig1_approximate_trace.csv) is a **visually approximated
trace of a published plot, NOT raw experimental observations or precise
experimental point digitization**. The approximately 5+10 ps explicit-water
pilot is **nonconverged**, does not match the source experiment's
protonation conditions, and does not validate a microscopic transition.
The pilot differs from the historical trace by hundreds of mV, even before
evaluating simulation convergence.

- [Full source/provenance and scientific interpretation](article/historical_benchmark.md)
- [Reproducible comparison generator](scripts/compare_historical_isotherm.py)
- [Integrity and no-fit regression tests](tests/test_historical_comparison.py)

Reproduce with `python scripts/compare_historical_isotherm.py` after
`pip install -e '.[plot,test]'`. The PDF workflow regenerates the figure
and PDF on source changes. A meaningful **predictive** comparison still
requires measured raw data, matched chemistry, converged multi-replica MD,
finite-size and electrostatic boundary validation, plus simultaneous
surface-pressure isotherms.
