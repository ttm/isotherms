# Molecular routes to Langmuir-film surface-potential isotherms: what is predictive, what is fitted, and what can be falsified?

**Research note (22 November 2026).** Companion to [`molecular_surface_potential.md`](molecular_surface_potential.md) and [`atomistic_surface_potential.md`](atomistic_surface_potential.md). This note compares several distinct molecular-modeling strategies; it is not an empirical validation. The run protocol in [`../scripts/explicit_water_pilot.py`](../scripts/explicit_water_pilot.py) is a short, neutral-acid, nonconverged experimental implementation. Any resulting voltage values are exploratory numerical observations, not evidence that the 30 Å² transition has been predicted.

## 1. Three distinct scientific targets

A complete theory must distinguish the **mechanical equation of state** $\pi(A,T,\mathrm{pH},I)$, the **electrical isotherm** $\Delta V(A,T,\mathrm{pH},I)$, and the **structural ensemble** $P(\text{chain tilt, conformation, hydration, charge, domain}\mid A,T,\mathrm{pH},I)$. It is possible to match any one without matching the others. The historical transition near 30 Å²/molecule is **not a physical constant**: temperature, composition, ionization, ions and compression history can shift phase boundaries and apparent onsets. For the specific historical data, the actual experimental conditions must be held fixed before comparison.

The comparison should be *forward*: select a molecular energy model from independent measurements or quantum calculations; calculate equilibrium structures and surface tension; reconstruct the potential from the **difference of complete neutral charge distributions** (film versus matched clean water) using the same boundary conditions. Do not select an area-dependent dielectric parameter, packing crossover, or headgroup bias to make $\Delta V(A)$ match a known curve.

## 2. Distinct modeling routes

| Route | Microscopic model | How $\Delta V(A)$ emerges | Strength | Main limitation | Priority |
|---|---|---|---|---|---|
| **A. Classical all-atom MD + explicit water** | Fixed-charge OPLS/CHARMM/OpenFF, neutral COOH and specified COO−/ions, explicit water | Planar charge-density difference + Poisson | Direct coordinates, structure, hydrogen bonds, pressure tensor; tractable | Force-field dependence, finite-size and sampling; fixed-charge induction absent | **1, initial benchmark** |
| **B. Polarizable MD / QM/MM** | Inducible dipoles or fluctuating charges, electronic headgroup QM | Self-consistent response and electrostatic profile | Tests whether fixed charges miss critical induction and water polarization | Large cost; electrostatic definitions and coupling need care | **2, mechanism check** |
| **C. Molecular density-functional theory (classical cDFT)** | Segment-resolved lipid density, orientational water density, ion densities, nonlocal interactions | Minimize molecular grand potential at fixed $A$, integrate resulting $\rho(z)$ | Analytic sensitivity, systematic thermodynamics, avoids long trajectories | Accurate nonlocal water and flexible-chain functionals challenging | **2, alternative predictive theory** |
| **D. Transferable statistical-mechanical lattice/cluster model** | Molecular packing and headgroup states with couplings extracted from MD/DFT, e.g. multistate Hamiltonian | Monte Carlo samples $P(s,w,\theta,\alpha)$; computed charge distributions | Cheap large domains, finite-size scaling, phase transitions | Calibration must be independent; an adjustable $A_0$ is not prediction | **1 in parallel to A** |
| **E. Molecular self-consistent field / polymer-chain theory** | Chain conformations, excluded volume, headgroup orientation, solvent dipoles | Solve Euler-Lagrange equations for conformation and charge profiles | Mechanism transparent, can scan phase diagram | Local correlation and H-bond structure may be oversimplified | **3** |
| **F. Martini / other coarse-grained MD + backmapping** | Fatty-acid and water beads, larger lateral patches | Emergent packing from coarse-grained interactions; recover atomic $\rho(z)$ only after backmapping and atomistic relaxation | Domain coexistence and hysteresis on large scales | Raw coarse-grained charges do not yield atom-resolved Kelvin potential; backmapping not unique | **3, large-scale structure** |
| **G. DFT/AIMD of local interfacial clusters** | Small headgroup-water-ion environments with electronic structure | Electron density + nuclei calculate local polarization; feed A/B/C/D | Headgroup protonation, induction, water response | Small-cell and short-time limits make direct large-scale isotherms difficult | **2, independent parameter anchors** |
| **H. ML interatomic potentials + molecular dynamics** | Local atomic energies trained to DFT for acid/water/ions | Sample configurations, coupled electrostatics from explicit electronic/charge model | Potentially longer simulations at higher accuracy | Long-range electrostatics, transferability, charging and surface physics must be validated | **research option** |

Routes are **not mutually exclusive**. An efficient program uses (A) for full-molecule observables, (G) to constrain uncertain interactions, (D) for sampling beyond small simulation cells, and (C) as a genuinely different thermodynamic theory. Route (B) tests whether predicted voltages survive a more realistic dielectric response. That combination is more informative than repeatedly refining curve fits.

## 3. One genuinely non-MD alternative: molecular density functional theory

Define an interfacial grand-potential functional $\Omega[\rho_f(\mathbf r,\omega),\rho_w(\mathbf r,\omega),\rho_\pm(\mathbf r)]$, where $\rho_f$ is the density of whole flexible fatty-acid conformers (orientation $\omega$ includes dihedrals), $\rho_w$ the orientational water density, and $\rho_\pm$ the ionic densities. A schematic construction is

$$
\Omega = F_{\mathrm{id}} + F_{\mathrm{packing}} + F_{\mathrm{chain}} + F_{\mathrm{HB}} +
F_{\mathrm{disp}} + F_{\mathrm{el}} - \sum_a\mu_a N_a,
\qquad
F_{\mathrm{el}}=\frac{\epsilon_0}{2}\int|\nabla\phi|^2\,d\mathbf r,
$$

with $-\epsilon_0\nabla^2\phi=\rho_{\mathrm{tot}}$ and the total (free + bound) molecular charge density built from site distributions rather than a fitted scalar interfacial permittivity. Charge regulation enters as a constrained reaction equilibrium, not as a post-hoc additive voltage. For a prescribed lipid coverage $\Gamma=1/A$, minimize $\Omega$ subject to fixed molecular count/area and electroneutrality. Compute $\Delta V(A)$ and $\pi(A)$ from the same optimum. Multiple local minima encode film phases; finite-area coexistence requires convexification/appropriate Maxwell construction. *Whether a particular functional correctly predicts the transition remains a test.* Incomplete hydrogen-bond and packing functionals can move the onset as surely as a poor force field.

## 4. A genuinely different small-molecule approach: cluster expansion from DFT

Calculate orientations and association energies of a small set of COOH/COO− groups, water bridges and neighboring chains with quantum chemistry under independently selected configurations. Derive site energies and multibody couplings, including headgroup rotation, hydrogen-bond occupancy and lateral adjacency. Sample a **two-dimensional statistical ensemble** in which the condensed fraction and characteristic area emerge from those interactions and geometry rather than from a set $A_0$. Couple it to a continuum ionic reservoir through an electrostatic free energy with neutralizing counterions. Its advantage is equilibrium and lateral-domain physics; its failure mode is hidden fitting in coarse-graining and inadequate representation of long-range Coulomb correlations.

## 5. An actual explicit-water pilot: meaning of the computation

The repository now includes `scripts/explicit_water_pilot.py`. It assembles two mirrored **neutral stearic-acid** monolayers around a TIP3P slab, uses OpenFF small-molecule parameters and OpenMM NVT dynamics, computes matched clean-water slab trajectories at each fixed $A$, and applies the existing planar-Poisson solver. Areas are supplied independently of any voltage data. The nominal grid is $24,28,30,32,36,45$ Å²/molecule. **No transition parameter exists in the builder or electrostatic operator.**

This is an integration/pilot and **not** sufficient to infer a transition from data because it has only 16 acids per interface, typically 5 ps of equilibration and 10 ps production, one seed, no verified pressure-area equation of state, fixed protonation, and a fixed-charge model. The water slab initially has lattice positions, *not an equilibrated water surface*. The structure must first relax fully. The two monolayer interfaces should agree within block/replica uncertainties; any mismatch is a red flag. Atomistic PME with vacuum and finite slab requires slab-size/field-correction checks, and force-field behavior in interfacial geometries requires validation. The **initial dataset should not be published as a predicted 30 Å² onset**, whatever curve it draws.

The coordinate-and-charge NPZ outputs use the same contract as the prior `atomistic_forward.py`. Inspection commands:

```bash
PYTHONPATH=src python scripts/explicit_water_pilot.py --areas 24 28 30 32 36 45 \
  --side 4 --equil-ps 5 --prod-ps 10 --interval-ps 1 \
  --out results/explicit_water_pilot

# Read each sample's algorithmically computed voltage (nonconverged):
cat results/explicit_water_pilot/summary.csv
```

To study the true molecular mechanism, repeat at **at least** lateral sizes 4×4, 8×8 and 12×12, independent seeds, multinanosecond equilibration and production, and at least two independently validated force-field/water combinations. Study approximately 10 mV level electrical convergence (provisional design target), plateau fields, slab vacuum, and finite-size scaling. The required simulation duration is decided by measured integrated autocorrelation times and slow domain reorganization, not an arbitrary fixed threshold; multiple nanoseconds may still be insufficient.

## 6. Experimental and numerical discrimination

1. **Define observables and sign:** $\Delta V_{\mathrm{film-clean}}$ measured at matched salt, pH, temperature and electrode convention. Compute matching periodic/slab electrostatic potential steps, not an absolute water potential. Confirm the neutral cell and integration origin independence.
2. **Match mechanical behavior:** derive surface tension from the pressure tensor or free-energy derivative with correct two-interface geometry; a voltage-only prediction is insufficient.
3. **Track structure:** condensed-domain fraction, chain tilt distribution, trans/gauche ratio, headgroup carbonyl tilt, hydrogen-bond network and first-shell water dipoles. Compare with GIXD, IRRAS, SFG, microscopy and any suitable independently measured scattering observables.
4. **Perturb chemistry:** change pH/salt with consistent charge regulation and explicitly compensating ions; then headgroup substitutions, chain lengths, temperature and compression speed. Treat hysteresis as a kinetic question.
5. **Compare model families:** A/B/C/D should predict *the same held-out structural and electrical observables*. If only one family can produce the transition with a tuned onset, the evidence for a molecular mechanism remains weak.
6. **Report null results:** a smooth $1/A$ potential, missing 30 Å² onset, reversed sign, or force-field-dependent transition is a scientifically useful falsification. Do not tune the crossover after seeing the curve.

## 7. Relevant literature (not claims of validation)

- Chachaj-Brekiesz, A. et al. (2022), *Comprehensive Approach to the Interpretation of the Electrical Properties of Film-Forming Molecules*, *J. Phys. Chem. B* **126**, 7037–7046, [doi:10.1021/acs.jpcb.2c04526](https://doi.org/10.1021/acs.jpcb.2c04526). DFT orientations and MD-informed three-layer regression, not an atomistic forward predictor of the stearic-acid isotherm.
- de Moura, A. F.; Oliveira Jr., O. N.; Trsic, M. (2005), *The Origins of the Surface Potential of a Stearic Acid Langmuir Film Investigated by Means of Molecular Dynamics Simulations*. Conference abstract, [SBQT proceedings](https://marcia.iqm.unicamp.br/SBQT13.pdf), earlier explicit-solvent stearic MD.
- *Structure and dynamics of water at liquid–vapour interfaces covered by surfactant monolayers of neutral stearic acid and charged stearate ions* (2008), [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0167732208000032). Water-interface comparison of different protonation states.
- *Molecular dynamic studies on Langmuir monolayers of stearic acid* (2013), [doi:10.1007/s40242-013-2301-x](https://doi.org/10.1007/s40242-013-2301-x). Simulations of the compression isotherm and acid–water hydrogen bonding.
- Gessner, C. (2018), *Stearic Acid Pressure-Area Isotherms*, Chapter 3 of doctoral thesis, [University of Reading](https://centaur.reading.ac.uk/78955/1/18009597_Gessner_thesis.pdf). Explicit CHARMM36 calibration and coarse-grained Martini phase behavior.
- *The analysis of neutron reflectivity from Langmuir monolayers of lipids using molecular dynamics simulations: the role of lipid area* (2025), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11919530/). Independent compression/expansion simulation strategies and structural constraints for a different lipid.
- OpenMM documentation, [running simulations](https://docs.openmm.org/latest/userguide/application/02_running_sims.html); [OpenMMForceFields](https://github.com/openmm/openmmforcefields). Software implementation descriptions, not evidence of correct interfacial physics.

## 8. Recommendation

Treat **molecular DFT + transferable lattice statistics** as the most promising *theoretical alternative* to atomistic MD, not as a replacement for independent structural data. Use all-atom MD as the immediate computational benchmark, *not* the only possible model. The decisive evidence will be whether independently parameterized A/B/C/D produce mutually compatible $\Delta V(A)$, $\pi(A)$ and structure without inserting the observed crossover by hand.