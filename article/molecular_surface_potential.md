# A coupled molecular-state and hydration model for the surface-potential isotherms of fatty-acid monolayers

**Status:** Theoretical research draft and computational demonstration (4 November 2026). No experimental parameters have yet been inferred from original measurements; no external experimental validation is claimed.

## Abstract

The surface-potential isotherm of a stearic-acid monolayer exhibits a low-potential expanded regime and a sharp increase upon compression to a molecular area of roughly 30 Å². Oliveira, Taylor and Morgan (1992) represented this change within a three-layer Demchak-Fort description by imposing a large reduction in an effective headgroup-region dielectric constant. Although phenomenologically effective, a prescribed permittivity discontinuity does not identify the microscopic degrees of freedom responsible for the measured voltage. Here we develop a coarse-grained molecular statistical model in which each fatty-acid molecule has an expanded or condensed packing state, a thermally distributed chain tilt, a carboxyl-group dipole and an associated interfacial-water orientational state. Intermolecular packing interactions generate correlated domains, while packing-dependent hydrogen-bond orientations determine excess normal water polarization. A separately and self-consistently charge-regulated Gouy-Chapman term represents dissociated acid groups and mobile subphase ions. The observable is the change in the normal dipole density relative to clean water plus the diffuse-layer contribution. No area-dependent *headgroup* dielectric parameter is introduced. A reproducible illustrative parameter set yields approximately zero potential at large molecular areas and a positive potential of 0.27 V at 20 Å², with a steep crossover near 30 Å². These results demonstrate mechanism feasibility, not a fit or proof of mechanism. We specify molecular observables and perturbations that could distinguish hydration-driven responses from tail-tilt, headgroup-conformation and electrostatic alternatives.

**Keywords:** Langmuir monolayers; stearic acid; surface potential; molecular statistical mechanics; hydrogen-bonded water; chain ordering; charge regulation; Gouy-Chapman; Demchak-Fort; molecular simulations.

## 1. Experimental motivation and scope

A monolayer surface potential is an interfacial electrical observable, conventionally the difference between the response of a covered aqueous surface and the corresponding clean subphase under the same conditions. A dipole interpretation must therefore account for the change in water polarization as well as the dipoles of the organic molecules. The measured potential must not be confused with the unmeasurable absolute microscopic electrostatic potential assigned to an arbitrary water model.

Three historical studies motivate our construction:

1. Taylor, Oliveira and Morgan (1989) found that a monolayer's carboxyl-group dissociation depends strongly on the concentration of monovalent salt. Their data supported Gouy-Chapman electrostatics *provided that charge regulation was taken into account*, with an apparent condensed-monolayer pKa of 5.7 +/- 0.2. They reported approximately 0.002 dissociation in ultrapure water and 0.034 in 0.01 M salt (their measurement conditions).
2. Taylor, Oliveira and Morgan (1990) reviewed the Demchak-Fort three-layer capacitor description and its limitations. They emphasized that headgroup, subphase-water and hydrophobic-terminal contributions are not cleanly recoverable from a single potential without additional assumptions. For condensed stearic acid they used a headgroup dipole of about 0.99 D, a terminal-group contribution near 0.33 D, and effective local permittivities of about 6.4 and 2.8.
3. Oliveira, Taylor and Morgan (1992) modeled the compression isotherm by a sharp change from bulk-water-like to low effective headgroup permittivity at about 30 Å²/molecule, accompanied by an assumed reorganization of water and domain formation. Their central claim is a plausible **interpretation**, not direct molecular observation of a local dielectric discontinuity.

**Prior molecular work exists.** Taylor and Bayes (1994) analyzed dipolar interactions beyond the basic capacitor approach. A 2005 conference abstract by de Moura, Oliveira Jr. and Trsic used OPLS-AA/SPC molecular dynamics of stearic-acid films, reproduced sensible packing but obtained a potential reportedly six times the experimental value. It resolved separate group contributions and oscillatory interfacial water polarization. Thus the present contribution cannot be described as the first molecular study of the problem. What we add is an explicit *area-dependent, charge-regulated molecular statistical ensemble* and a program of joint electrical and structural discrimination. Whether that is publishably novel relative to the full subsequent literature remains to be established by a systematic review.

## 2. Microscopic degrees of freedom

Consider `N` C18 acid molecules in a surface area `S = N A`, with `A` the prescribed area per molecule. A molecular site `i` has:

- a packing state `s_i` in `{0,1}`, expanded or condensed;
- a continuous chain inclination `u_i = cos(theta_i)` in `[0,1]`, with theta measured from the surface normal;
- a local hydration-water orientation `w_i` in `{-1,+1}`, defining the sign of the **excess** water dipole along the normal relative to clean water;
- an acid dissociation state `d_i` in `{0,1}`, meaning COOH or COO-.

The `w_i` state is a *coarse-grained first-shell water orientational mode*, not a claim that bulk water has only two Euler orientations or that there is exactly one physical water molecule per acid. Its dipole amplitude absorbs the number and correlated orientation of water molecules whose configuration changes upon adsorbing an acid molecule. Hydration vacancy, many-water orientational continua and lateral hydrogen-bond bridges are natural extensions.

One economical dimensionless Hamiltonian for packing and hydration is

\[
\beta H_{\mathrm{mol}} =
  \delta(A)\sum_i s_i
  - J_c\sum_{\langle i,j\rangle}\sigma_i\sigma_j
  -\sum_i \bigl[B_e(A)+\lambda s_i\bigr]w_i
  - J_w\sum_{\langle i,j\rangle}w_i w_j
  + \sum_i \frac{\kappa_{s_i}}{2}(u_i-u^0_{s_i})^2,
\qquad \sigma_i=2s_i-1.
\]

Here `beta = 1/(k_B T)`. The field `delta(A) = g_A (A-A_0)` is the free-energy preference for condensation induced by compression. Its crossover area `A_0` is a **model parameter**, not a value calculated ab initio. The molecular packing coupling `J_c` favors neighboring molecules in the same phase. `B_e(A)` is the orientation-dependent free energy of associated water in an expanded local environment; `lambda` modifies it in a condensed environment. `J_w` represents lateral cooperativity in the hydration network. For the analytical mean-field curves we take `J_w = 0`, allowing exact tracing of the two water states conditional on packing. The supplied code also supports nonzero `J_w` in the Monte Carlo sampler.

The packing states and water orientations are **molecular variables**, not hidden switches in a local permittivity. The model does retain coarse-grained phenomenological free-energy parameters, so it is not an ab initio prediction.

## 3. Mean-field reduction and domain sampling

For `J_w=0`, summing over `w_i = +/-1` at fixed packing state produces the local water partition functions

\[
Z_w(s)=2\cosh[B_e(A)+\lambda s],\qquad
\delta_{\mathrm{eff}}(A)=\delta(A)-\ln\frac{\cosh[B_e(A)+\lambda]}{\cosh B_e(A)}.
\]

With `phi = <s_i>` and four neighbors per molecular site, the dimensionless mean-field free energy per molecule (up to a `phi`-independent term) is

\[
 f(\phi) = \phi\,\delta_{\mathrm{eff}}
 -2J_c(2\phi-1)^2
 +\phi\ln\phi+(1-\phi)\ln(1-\phi).
\]

The equilibrium condensed fraction minimizes `f(phi)` over `[0,1]`; this specification is more robust than choosing an arbitrary branch of the self-consistency equation. The corresponding excess hydration order parameter is

\[
\bar w=(1-\phi)\tanh B_e(A)+\phi\tanh[B_e(A)+\lambda].
\]

A Gaussian orientational energy on the bounded interval `0 <= u <= 1` yields

\[
\bar u_s=\frac{\int_0^1 u\exp[-\kappa_s(u-u_s^0)^2/2]\,du}{\int_0^1\exp[-\kappa_s(u-u_s^0)^2/2]\,du},\qquad
\bar u=(1-\phi)\bar u_e+\phi\bar u_c.
\]

This is thermally averaged molecular chain tilt, not geometrical inference from a conical volume. A separate two-dimensional checkerboard Metropolis sampler (`metropolis_point`) retains finite-size spatial correlations in the `s_i` and `w_i` fields. The deterministic mean-field curves and sampled lattice maps answer different questions. In particular, a rapid crossover is **not automatically a thermodynamic phase transition**. At the illustrative `J_c=0.20`, the square-lattice Ising coupling is below its zero-field critical coupling, and the mean-field energy also lacks a double-well at its symmetry point. Claiming a true first-order transition would be unjustified.

## 4. Charge regulation and diffuse-layer electrostatics

For a monovalent symmetric electrolyte of ionic strength `I` mol/L, let `n_b = 1000 N_A I` be its bulk number density in m^-3. The surface charge density is `sigma = - e alpha / A`, where `alpha = <d_i>` is the mean ionized fraction and `A` is converted to m². For a planar Gouy-Chapman layer,

\[
 \psi_0=\frac{2k_BT}{e}\,\operatorname{asinh}\!\left(
 \frac{\sigma}{\sqrt{8\epsilon_0\epsilon_w k_BT n_b}}
 \right).
\]

For the reaction COOH -> COO- + H+, the local proton concentration at a negative interface exceeds that in the bulk. Hence

\[
\ln\frac{\alpha}{1-\alpha}=\ln 10\,[\mathrm{pH}_{\mathrm{bulk}}-\mathrm{p}K_a^\mathrm{eff}]
 + \frac{e\psi_0}{k_BT}.
\]

These two equations are solved simultaneously at each area. The **plus sign** in front of `e psi_0` is essential: since `psi_0 < 0`, self-repulsion of ionized carboxyl groups *suppresses* dissociation. The analytic implementation allows optional changes of apparent pKa with packing and hydration, but defaults them to zero because such couplings have not been independently estimated. The current implementation solves structural order first and charge regulation second. When the optional pKa shifts are activated, this is a one-way coupling rather than a fully variational mutual feedback; extending the Hamiltonian by the diffuse-layer charging free energy and re-minimizing all order parameters is required for thermodynamic self-consistency.

`epsilon_w=78.5` is the **bulk** water relative permittivity in continuum diffuse-layer electrostatics. That use does not resurrect the hypothesis of a discontinuous *local headgroup* permittivity. At large ionic strength or large local electric field, finite ion size, Stern layers, non-Poisson-Boltzmann correlations and field-dependent dielectric response must be considered.

## 5. Surface potential from microscopic polarization

The fundamental electrostatic observable follows from the planar average of the excess charge density. For an electrically neutral complete simulation cell, with vanishing far-field electric fields, the potential difference across an interface is

\[
\Delta\varphi = \frac{1}{\epsilon_0}\int z\,\Delta\rho(z)\,dz.
\]

`Delta rho` must include all partial charges, water, headgroup atoms, chains and compensating ions, and must be evaluated using the **same electrostatic convention** for the monolayer and clean-water reference. An absolute water-surface potential computed from a particular force field or quadrupole partition cannot directly be equated to the measured Kelvin-probe change. Polarizable molecular models and slab boundary conditions require particular care. Charge-density integration is preferable for future atom-resolved validation.

In the current coarse-grained dipolar reduction the expected excess normal dipole per fatty-acid molecule is

\[
\begin{aligned}
\mu_\perp(A) ={}& (1-\alpha)\mu_{\mathrm{COOH}}+\alpha\mu_{\mathrm{COO^-}} \\
&+\mu_{\mathrm{CH_3}}\bar u(A)+\mu_w^{\mathrm{ex}}\bar w(A).
\end{aligned}
\]

The modeled observable is

\[
\boxed{\quad \Delta V(A) = \frac{\mu_\perp(A)}{\epsilon_0 A}+\psi_0[A,\alpha(A)].\quad}
\]

The quantity `mu_w^ex * wbar` represents the water polarization change *relative to clean water*, including molecules displaced or reoriented by acid adsorption. The plus sign for an upward dipole is conventional here, with the subphase-to-air normal positive. The decomposition is an interpretable **model partition**; it is not an experimentally unique separation. The historical formula instead defines screened contributions `mu_i/epsilon_i`, which are not interchangeable with unscreened local molecular dipole moments.

For `A=20 Å²`, even a single 1 D normal dipole corresponds to roughly **1.88 V** in the unscreened Helmholtz relation. Therefore cancellation by associated-water polarization and the ionic diffuse layer is large. A 0.2 D error in the molecular excess dipole can move the predicted voltage by approximately 0.38 V. Reproducing a 0.25 V signal alone does not uniquely determine water order, molecular polarizability or local permittivity. This is a central identifiability challenge.

## 6. Computed example and parameter provenance

The implementation includes a small **illustrative**, hand-selected parameter set. It was selected to test whether the mechanism can generate a qualitatively appropriate surface-potential isotherm, and must not be interpreted as independently known physical constants or a best fit to the 1992 points.

| Parameter | Illustrative value | Interpretation |
|---|---:|---|
| `T` | 295 K | experimental reference temperature |
| `pH`, `pKa` | 5.6, 5.7 | based on the 1989 study's conditions/estimate |
| `I` | 2.5 x 10^-6 M | effective dilute monovalent electrolyte condition |
| `A_0` | 30 Å² | packing crossover set near historical onset, **not predicted** |
| `g_A` | 0.45 (kBT)/Å² | packing free-energy slope, **assumed** |
| `J_c`, `J_w` | 0.20, 0 | dimensionless packing and water-neighbor couplings |
| `B_e(A)` | `-0.37 - 0.12(40/A - 0.5)` | **assumed** area-dependent water orientation field |
| `lambda` | +0.09 | change in hydration field for condensed packing |
| `mu_COOH`, `mu_CH3` | 0.99 D, 0.33 D | historical group-dipole estimates, convention dependent |
| `mu_COO-` | 0.75 D | **illustrative**; negligible at low ionization |
| `mu_w^ex` | 2.5 D | **illustrative effective** hydration-water dipole amplitude |
| `u_e^0`, `u_c^0` | 0.48, 1.02 | centers of bounded thermal tilt distributions; `u_c^0 > 1` centers a truncated distribution near upright |
| `kappa_e`, `kappa_c` | 8, 150 | angular stiffness, dimensionless and assumed |

Selected deterministic mean-field outputs (rounded, computed by `scripts/run_model.py`):

| Area (Å²) | Condensed fraction | Excess water order | `psi_0` (mV) | Total `Delta V` (mV) |
|---:|---:|---:|---:|---:|
| 80 | 0.000 | -0.354 | -126 | -1 |
| 40 | 0.002 | -0.405 | -138 | -8 |
| 35 | 0.022 | -0.418 | -140 | -23 |
| 30 | 0.455 | -0.404 | -143 | 120 |
| 28 | 0.893 | -0.381 | -144 | 303 |
| 25 | 0.977 | -0.392 | -146 | 332 |
| 20 | 0.998 | -0.430 | -150 | 273 |

The sharp increase results from **coupled changes in chain orientation and hydration-water polarization when condensed packing becomes prevalent**, plus the geometrical `1/A` factor. The model also produces a shallow condensed-regime maximum. None of these numbers is a digitized experimental observation. The shape is not sufficient evidence for the water mechanism since the water bias and crossover area were chosen illustratively.

## 7. How this changes Demchak-Fort

The earlier three-layer model is

\[
\Delta V_{DF}=\frac{1}{\epsilon_0 A}
\left(\frac{\mu_1}{\epsilon_1}+\frac{\mu_2}{\epsilon_2}+\frac{\mu_3}{\epsilon_3}\right)+\psi_0.
\]

Our model retains the physical distinction between water, head and tail, but replaces its area-dependent effective local dielectric factors with **ensembles of molecular orientations** and explicit charge regulation. This does not prove the existence or nonexistence of a low interfacial dielectric susceptibility. Indeed, the local differential dielectric response of orientationally constrained water can be low even if the average hydration polarization is large. In a polarizable extension one would compute dielectric response from polarization fluctuations, for example a properly boundary-corrected normal susceptibility, rather than deducing it from a fit to a voltage curve. A single scalar `epsilon_2` is not guaranteed to be the appropriate descriptor of an anisotropic, molecularly thin interface.

The historical reconstruction in `src/isotherms/legacy.py` uses the published approximate contributions and a literal step in epsilon2 at 30 Å². It is an implementation of the model assumptions, not a fit to digitized measurements; differences near the discontinuity also reflect idealized threshold handling.

## 8. Discriminating experimental tests

**Joint measurements are essential.** A voltage curve alone is underdetermined. The following predictions or tests can break some of the degeneracy:

1. **Salt/pH perturbation at fixed packing:** as monovalent ionic strength increases, the model's negative diffuse-layer potential decreases in magnitude while the ionization fraction increases. Measure `Delta V`, surface pressure and independently estimated charge simultaneously. Changes of packing with salt would require extending the structural Hamiltonian before attributing all variation to Gouy-Chapman.
2. **Structure/voltage co-registration:** use Brewster-angle microscopy or grazing-incidence X-ray diffraction for condensed area fraction and structure, infrared reflection-absorption spectroscopy for chain tilt, and surface potential at the same area and compression rate. The model predicts that the steep electrical change approximately follows the increase of `phi` and `ubar`. A potential jump without structural change would challenge the assumed mechanism.
3. **Hydration orientation:** use interface-selective vibrational sum-frequency generation or comparable techniques sensitive to hydrogen-bonded water orientation. Compare perturbations of water orientation with `Delta V` after correcting for measured chain tilt and surface charge. The two-state water approximation predicts a correlated sign and amplitude, but optical response is not directly proportional to `wbar` without spectral modeling.
4. **Water-network disruption:** alter subphase ion identity or temperature while independently checking ionic strength and pH. The 1992 work reports that SCN- changes the onset area. A mechanism involving hydration networks predicts sensitivity beyond diffuse-layer screening, but neither the **direction nor the size** of the shift follows uniquely from this preliminary Hamiltonian; ion-specific adsorption must also be considered.
5. **Chain-length and terminal substitutions:** varying the alkyl chain or distal group alters `mu_CH3`, orientational stiffness, dispersive packing and possibly the hydration structure. Independently measure tilt and charge to determine whether a transferable water parameterization is possible.
6. **Compression/expansion rate and hysteresis:** spatially correlated Monte Carlo ensembles can be extended to dynamics, nucleation barriers, defects and compression protocols. A quasi-equilibrium prediction should not be compared to an unrelaxed experimental isotherm without a kinetics model.

A strong falsification criterion is a parameter set calibrated on a limited number of conditions that fails to predict **held-out** pH, ionic-strength, temperature or chain-length trends, or that requires mutually inconsistent hydration orientations across measurements.

## 9. Parameter estimation, reproducibility and uncertainty

The Python implementation uses immutable parameter records, an exact single-site trace for the two hydration orientations when `J_w=0`, free-energy minimization for condensed fraction, bounded quadrature for tilt, a root solve for ionization, and a periodic two-dimensional checkerboard Metropolis sampler. All figure scripts print a warning that the curves are illustrative. The test suite checks physical bounds, charge-regulation consistency, salt-screening response, Monte Carlo repeatability, and benchmark behavior.

The optional `scripts/fit_data.py` accepts original data with `area_a2,potential_mv,sigma_mv` and estimates only three parameters (`A_0`, the expanded-phase water bias, and condensed-phase hydration shift), while holding others fixed. It records chi-square and a *local* least-squares covariance. Such a covariance does not include uncertain dipole assignments, systematic voltage zero shifts, correlations along a compression cycle, calibration drift or structural model inadequacy. Preregistered train/test splits and multiple independent datasets are preferable to optimizing all free energies simultaneously. No raw measurements have been digitized or invented in this repository.

Further verification should compare the mean-field predictions against larger lattice Monte Carlo simulations over multiple random seeds, finite-size scaling, autocorrelation times, and equilibrium histories. It should also compare the baseline with ablations (`mu_w^ex=0`, fixed tilt, no ionization, `J_c=0`) to isolate which microstates generate which part of the voltage.

## 10. Physical and conceptual limitations

The current theory is intentionally a **minimal molecular-state hypothesis**, not a closed quantitative molecular solution. Its main limitations are:

- The packing preference `delta(A)` and crossover area are prescribed. They do not yet emerge from a molecular force field or a computed pressure-area equation of state.
- A single excess-water orientation spin per headgroup is too small a state space to resolve multilayer hydrogen bonding, angular correlations, water quadrupoles, displacement, hydration occupancy or proton hopping.
- The headgroup's vertical dipole is fixed except for protonation, whereas carbonyl rotations and local acid conformers are possible alternative explanations of the voltage anomaly.
- Electrostatic polarization of headgroups and water by neighbors is only represented indirectly in effective moments and local water energies; the Coulomb and induction energies are not solved explicitly.
- The diffuse layer is a mean-field planar Gouy-Chapman continuum; at sufficiently high ionic strength or charge density it is incomplete.
- The molecular and ionic sectors are currently solved sequentially, not through a fully coupled variational Gibbs free energy when packing-dependent pKa corrections are enabled.
- A macroscopic vibrating-plate/Kelvin measurement and an atom-resolved electrostatic calculation must be matched by their **potential difference relative to clean water**, with compatible conventions and treatment of potential offsets.
- The model is underidentified by a single potential isotherm; many distinct microscopic scenarios can yield similar voltage-area curves.

## 11. Road map to a fully atom-resolved test

A rigorous next stage would use explicit-water molecular dynamics with stearic acid and stearate at experimentally relevant areas and carefully controlled ionization, ideally polarizable water/headgroup force fields or ab initio benchmarks for local dipole response. For each frame: (i) record all site and ion coordinates and partial charges; (ii) compute planar charge-density differences against a clean-water control with the same model, boundary convention, slab electrostatics and thermodynamic conditions; (iii) integrate the planar Poisson equation rather than infer local epsilon from a fitted voltage; (iv) separately analyze carbonyl orientations, chain tilt, hydration-shell occupancy, water orientation, hydrogen bonds and aggregate domains; (v) assess polarization convergence, lateral size, ionization, simulation time and force-field dependence. Coarse-grain the resulting structural distributions back to the lattice model, then see whether fitted parameters transfer across pH, salt and compression state. This would directly probe whether interfacial hydration reorganizes near the historical onset and could also reject the present model.

## 12. Conclusion

A microscopic mechanism exists that can qualitatively generate the long-observed stearic-acid potential-area response without *postulating* a step in a local headgroup dielectric constant. Compression changes the relative populations and correlations of molecular packing states, chain orientations and interfacial-water orientations; their excess polarization is converted to voltage while the degree of acid ionization responds to the diffuse-layer potential. The immediate scientific result is a well-defined, runnable and falsifiable **candidate mechanism**. It is not evidence that hydration reorientation is dominant, not an independent prediction of the 30 Å² crossover, and not yet a validated replacement for Demchak-Fort. Competing microscopic mechanisms, independent structural measurements and genuine held-out isotherms must determine its scientific status.

## References

1. Taylor, D. M.; Oliveira Jr., O. N.; Morgan, H. *The surface potential of monolayers formed on weak acidic electrolytes: implications for lateral conduction*. Chemical Physics Letters **161** (1989), 147-150.
2. Taylor, D. M.; Oliveira Jr., O. N.; Morgan, H. *Models for interpreting surface potential measurements and their application to phospholipid monolayers*. Journal of Colloid and Interface Science **139** (1990), 508-518.
3. Oliveira Jr., O. N.; Taylor, D. M.; Morgan, H. *Modelling the surface potential-area dependence of a stearic acid monolayer*. Thin Solid Films **210/211** (1992), 76-78. https://doi.org/10.1016/0040-6090(92)90172-8
4. Demchak, R. J.; Fort Jr., T. J. *Surface dipole moments of close-packed un-ionized monolayers at the air-water interface*. J. Colloid Interface Sci. **46** (1974), 191-202. https://doi.org/10.1016/0021-9797(74)90002-2
5. Oliveira Jr., O. N.; Taylor, D. M.; Lewis, T. J.; Salvagno, S.; Stirling, C. J. M. *Estimation of group dipole moments from surface potential measurements on Langmuir monolayers*. J. Chem. Soc. Faraday Trans. 1 **85** (1989), 1009-1018. https://doi.org/10.1039/F19898501009
6. Taylor, D. M.; Bayes, G. F. *Calculating the surface potential of unionized monolayers*. Physical Review E **49** (1994), 1439. https://doi.org/10.1103/PhysRevE.49.1439
7. de Moura, A. F.; Oliveira Jr., O. N.; Trsic, M. *The Origins of the Surface Potential of a Stearic Acid Langmuir Film Investigated by Means of Molecular Dynamics Simulations*. Abstract, XIII SBQT (2005), https://marcia.iqm.unicamp.br/SBQT13.pdf (conference proceedings, not a full peer-reviewed article).
8. *Molecular Dynamics Simulation of Hydrated DPPC Monolayers Using Charge Equilibration Force Fields*. 2012. https://pmc.ncbi.nlm.nih.gov/articles/PMC3488352/ (example of planar charge-density integration, not validation of stearic-acid predictions).
