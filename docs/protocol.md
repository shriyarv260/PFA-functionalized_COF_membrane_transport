# Scientific protocol and decision gates

## 1. Define a controlled design

Start with one stable, experimentally supported COF chemistry, one stacking arrangement and a matched parent/functionalized pair. Expand to at least three accessible pore sizes only if chemically distinct, stable structures exist. Proposed surface groups include neutral hydroxyl and amine groups and a positively charged site; the latter requires a specific protonation model and counterions. Do not add abstract charges without a molecular definition. Evaluate seven PFAS species separately before mixtures; use at least three independent seeds per condition. Record temperature, ionic strength, solute concentration, pH assumption, membrane thickness and accessible area.

## 2. Build and validate systems

Use CURATED-COF IDs to trace original crystallographic references; prefer optimized structures from Materials Cloud after inspecting both. Check geometry, bond orders, disorder, connectivity, stacking and aqueous stability. Construct a periodic slab with water reservoirs along z; terminate severed bonds chemically. Compare rigid and flexible treatments if framework motion is constrained. Quantify pore-limiting diameter on the relaxed structure using a declared probe radius and analysis method.

Parameterize the framework, PFAS anions, water and ions using one compatible nonbonded convention. Document charges, torsions, fluorocarbon cross-interactions and combining rules. A CIF or a PubChem SMILES is not a force field. Validate against available hydration, bulk diffusion, density, structural and adsorption evidence before using a generic force-field assignment. Neutral-acid database properties cannot substitute for anion charges/topologies.

## 3. Equilibration and unbiased MD

Minimize, equilibrate solvent with restrained framework, and then equilibrate the chosen membrane ensemble. Select pressure coupling according to slab geometry; do not apply isotropic volume coupling to an unsupported membrane/vacuum setup. Inspect density, energy, temperature, pressure, membrane spacing and PFAS position. Use independent starts and discard equilibration based on observables. The provided 20 ns production duration is an initial planning value, not a convergence guarantee.

For MSD, make molecules whole, unwrap their trajectories, and remove global translation using the membrane reference. Record any rotational alignment. Fit in the diffusive regime and vary the fit interval. Report in-plane diffusion separately; a confined normal-coordinate MSD can plateau. Estimate uncertainty across independent replicas, not by ordinary least-squares errors over correlated lags. A local D(z) requires a validated restrained-dynamics estimator and is not implemented here.

## 4. Umbrella sampling

Define signed z displacement between a membrane reference group and one PFAS COM. Span bulk water on both sides, the interfaces and pore interior; select any lateral pore restraint explicitly and include its free-energy consequences. Use starting configurations with the molecule physically in each window. Equilibrate each umbrella before collecting statistics; retain a common coordinate definition and restraint convention.

Initial planning values: 0.1 nm spacing and 1000 kJ/mol/nm² springs. Tune them using actual overlap and relaxation; neither is universally suitable. Choose each pull group's PBC reference atom and audit wrapping. The box must support the full coordinate range without image ambiguity. Check orthogonal relaxation, hydration transitions, pore occupancy, hysteresis, window histograms and autocorrelation. Repeat sampling from both directions or independently seeded starts. Reconstruct PMFs using consistent temperature, bins and reference bulk region; quantify bootstrap and replica uncertainty. Numerical WHAM convergence alone is insufficient.

```bash
gmx wham -it tpr-files.dat -ix pullx-files.dat -o pmf.xvg -hist histograms.xvg -unit kJ -temp 298.15 -b 5000 -ac -nBootstrap 200 -bsres bootstrap.xvg
```

The 5000 ps discard is a starting example requiring equilibration evidence. Restrict input lists to one common physical system; never combine different species or membrane chemistries into one PMF.

## 5. Additional observables

Hydration: compute headgroup–water oxygen RDFs and coordination numbers. Select cutoff from an RDF minimum and inspect sensitivity; the helper uses an explicit orthorhombic minimum-image calculation. Residence: define an adsorption region or contact criterion, compare continuous and intermittent definitions, and use censor-aware survival analysis. The event extractor flags censoring; it does not supply a fitted lifetime. Molecules already adsorbed at the first frame are left-censored.

Water: count complete reservoir-to-reservoir passages using molecule identities and two boundaries. Estimate steady net flux across multiple pressure differences, with the effective hydraulic-minus-osmotic driving pressure and a declared molecular volume. Report permeance in L m^-2 h^-1 bar^-1. Multiply by membrane thickness only when explicitly reporting thickness-normalized permeability. Do not equate equilibrium bidirectional crossings with pressure-driven flux.

## 6. Structure–property conclusions

Join species, actual pore diameter, group chemistry/density, charge, hydration, PMF barrier relative to bulk, adsorption well depth, transport and uncertainty by system/replica IDs. Use matched contrasts and confidence intervals; avoid overfitting a small matrix. Report observed rejection only from steady-state feed/permeate concentrations with concentration polarization and adsorption transients addressed. PMF-based relative transport requires a model, standard-state/area conventions and D(z); barrier-only ratios are not measured rejection.

## Completion criteria

A publishable result needs validated force fields, deposited structures, sufficient independent sampling, connected overlapping windows, replica agreement, stable fit intervals, and an auditable link from raw trajectories to every plotted value. The current repository establishes the analysis foundation; none of these production gates is claimed complete.
