# Explicit-water stearic-acid pilot: numerical and scientific audit (22 November 2026)

**Reproducibility status:** This audits the completed GitHub Actions [MD pilot run](https://github.com/ttm/isotherms/actions/runs/45057869556) and its trajectory artifact (`neutral-stearic-acid-explicit-water-pilot`, artifact ID 14769025865), together with the committed `results/explicit_water_pilot/summary.csv` and `protocol.json`. The run passed as software, but the output does **not** pass scientific convergence criteria. Numerical quantities below are measured on the actual archived trajectory files, not invented points or experimental observations.

## 1. Reproducibility and CI

- Main-line tests at [run 45057869567](https://github.com/ttm/isotherms/actions/runs/45057869567): **15 passed, 1 skipped**. The skipped module is the optional RDKit geometry test in the pip-only environment. The MD workflow executes that module separately: **8 passed**.
- Explicit-water workflow [run 45057869556](https://github.com/ttm/isotherms/actions/runs/45057869556): **success**, six area-specific runs; results published by bot commit. Earlier runs failed due to an OpenMM API incompatibility and were corrected; historical failures remain in Actions.
- The setup includes 16 neutral C18H36O2 molecules per interface, two interfaces, a 4 nm initially lattice-filled TIP3P water slab, 16 nm box height, OpenFF 2.2.1, 3D PME, 295 K, fixed-area NVT, 5 ps nominal equilibration, 10 ps nominal production, ten snapshots at 1 ps spacing, and one initialized configuration per area. **There is no explicit variable protonation, pressure-area calculation, independent replicate or calibrated electrical reference.**
- The CI success demonstrates that the molecular builder, force-field assignment, OpenMM integration, Poisson postprocessor, and archive/commit pipeline can execute on the runner. It does not establish thermodynamic equilibrium or correct physical voltages.

## 2. Numerical audit of archived trajectories

For each independently generated film and clean-water trajectory, the periodic planar-Poisson operator was rerun **on each stored frame**, with the same 512 bins and water/air windows used by the workflow. `Frame SD` is the sample standard deviation of the *10 correlated, paired-by-index film-minus-clean snapshots* and is **not a confidence interval, independent standard error, or block-averaged uncertainty**. `Side difference` uses the two interfaces' *time-averaged* potential differences. `Water-window span` is the larger of the film and clean potential ranges inside the nominal 7.6-8.4 nm bulk-water sampling interval, computed from their time-averaged profiles.

| Area (Å²/molecule) | Mean ΔV (V) | Frame SD (V) | Side difference (mV) | Max water-window span (mV) |
|---:|---:|---:|---:|---:|
| 24 | 1.480 | 0.555 | 122.9 | 255.1 |
| 28 | 1.123 | 0.403 | 116.9 | 272.4 |
| 30 | 0.917 | 0.517 | 51.1 | 232.3 |
| 32 | 1.090 | 0.563 | 54.4 | 209.2 |
| 36 | 0.717 | 0.180 | 79.2 | 183.5 |
| 45 | 0.469 | 0.491 | 42.6 | 166.3 |

**Interpretation:** All six nominal bulk-water windows are non-flat by over 0.10 V in at least one system. The frame SDs (0.18-0.56 V), side asymmetries (43-123 mV), and nonmonotonic values at 30 versus 32 Å² are too large to identify an equilibrium onset at 30 Å² or to claim an agreement with historical surface-potential isotherms. The fact that all frames remain charge-neutral is necessary but not sufficient to validate their electrostatics. Real bulk-water and adjacent-air plateaus, verified slab corrections, and autocorrelation-aware errors are prerequisites. The results may be artifacts of unequilibrated water and headgroups, very small lateral geometry, limited statistics, and force-field/PME limitations. These explanations are plausible, **not proven by this pilot**.

## 3. Scientific and implementation quality gates

| Gate | Present status | Minimum credible next test |
|---|---|---|
| Simulation launches and outputs correct shapes | **Pass** | Preserve in CI with exact package and platform versions |
| Charge neutrality and analytic Poisson-sheet benchmark | **Pass for tested cases** | Add realistic slab and box-size numerical benchmarks |
| Correct electrostatic plateau extraction | **Not established** | Verify flat, stable water and air regions for every area, frame block and replica |
| Equilibrium and sampling convergence | **Fail / not demonstrated** | Multinanosecond-or-longer trajectories guided by autocorrelation, independent seeds and equilibration checks |
| Thermodynamic mechanical isotherm π(A) | **Absent** | Surface tension/pressure calculation or free-energy derivative, with two-interface conventions |
| Microscopic interpretation | **Absent** | Chain tilt, headgroup dihedrals, hydrogen bonds, water orientation/occupancy, and domain statistics |
| Physical transferability and force-field robustness | **Absent** | Independent force field/water comparison, slab correction and lateral-size scaling |
| Protonation/pH/salt response | **Absent** | Neutral benchmark first, then explicit correctly neutralized stearate/ion ensembles |
| Experimental validation | **Absent** | Original measured data with pH, salt, T, compression rate, sign and uncertainty, held out from parameter selection |

## 4. Documentation, code, and paper review

The original 2026 molecular-state manuscript explicitly selects a packing crossover parameter near 30 Å²; it proves *mechanistic possibility*, not ab initio predictive onset. The 22 November atomistic paper correctly defines a molecular forward problem and provides a tested Poisson solver, but it was initially written before the short MD pilot was run. Its text has been updated to distinguish an actual **nonconverged pilot** from missing **converged production trajectories**, and its PDF benchmark caption has been deduplicated. The original isotherm manuscript has not been altered.

The `src/isotherms/atomistic.py` operator is useful for neutral periodic slabs and analytic validation. However, the library-level `excess_surface_potential` comparison only inspects the z-grid, not the lateral box dimensions; the higher-level CLI checks all three box lengths. This is a **library API safeguard gap** that should be closed when refining the analysis data contract. The MD tool does not yet calculate validation metrics or stop on non-flat sampling windows. Package versions, initial structure, electrostatic corrections, and hardware should be recorded explicitly for archival reproduction. The main README previously described the atomistic calculation as a proposal; this audit supersedes that historical status statement.

**Consequence:** Current `main` is an executable methods/prototype repository, suitable for research development and numerical benchmarking. It is **not** a validated molecular prediction of the experimental surface-potential isotherm and should not be represented as one.