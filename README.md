# Atomistic Modeling of PFAS Transport in Functionalized COF Membranes

Reproducible public-data assessment, molecular-simulation planning, and transport analysis for PFAS–COF research.

**Status:** working analysis software and public molecular/COF metadata. No original trajectories were provided and no production atomistic simulations have been run. Synthetic verification outputs are explicitly labeled. This repository was assembled in September 2026; the supplied May–August 2025 project dates are historical context, not dates of these calculations.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install .
python -m unittest discover -s tests -v
pfas-cof demo --out results/synthetic
python scripts/analyze_public_data.py
```

The included public snapshot supports offline analysis. To refresh public records (network required):

```bash
python scripts/fetch_public_data.py
python scripts/analyze_public_data.py
```

The COF source revision is pinned in `data/sources.lock.json`; PubChem response bytes and hashes are retained because live properties can change.

## Research questions

How do accessible pore dimensions, hydration and chemically explicit functional groups affect PFAS partitioning, residence and transport? How do carboxylate/sulfonate headgroups and chain length change those responses? What are the tradeoffs between water permeance, PFAS passage and adsorption?

The proposed panel is PFBA, PFHxA, PFOA, PFNA, PFBS, PFHxS and PFOS. Public properties refer to parent acids; proposed aqueous simulations use explicitly parameterized anions under a documented solution condition. Pore sizes and surface chemistries must be measured on relaxed structures, not inferred from names.

## What is included

| Path | Content |
|---|---|
| `data/public/` | Real PubChem responses and CURATED-COF registry, upstream license, checksums |
| `results/public/assessment.md` | Descriptive assessment and missing-evidence inventory |
| `src/pfas_cof/` | MSD/diffusion, harmonic WHAM, hydration, censored residence episodes, water permeance |
| `scripts/` | Data acquisition, public-data report, umbrella input generation |
| `simulation/` | GROMACS base settings and workflow requiring validated atomistic inputs |
| `docs/` | Scientific protocol, data contract, sources and interpretation limits |
| `tests/` | Mathematical and boundary-condition validation |
| `results/synthetic/` | Seeded estimator verification, not PFAS transport predictions |

## Analyze real data

Export aligned, unwrapped molecular COM coordinates in nm to an NPZ containing `positions_nm` (frames × molecules × 3) and scalar `dt_ps`. Use only equilibrated, unbiased trajectory segments for MSD.

```bash
pfas-cof msd trajectory.npz --axes xy --max-lag 1000 --fit 100 500 --out results/real/diffusion
pfas-cof wham windows.csv --range -3 3 --bins 120 --temperature 298.15 --out results/real/pmf
```

The umbrella manifest uses `path,center_nm,k_kJ_mol_nm2`; each path points to a one-column coordinate file in nm, relative to the manifest. Remove equilibration and subsample according to measured autocorrelation first. The reference WHAM implementation checks adjacent histogram connectivity, preserves empty bins as missing, and fails on nonconvergence. For production use, compare to GROMACS `gmx wham`, estimate autocorrelation, bootstrap, and check independent replicas.

See [protocol](docs/protocol.md), [data contract](docs/data_contract.md), and [sources](docs/sources.md). A high PMF barrier alone does not establish filtration rejection. Confinement may invalidate a global Einstein diffusion fit; use in-plane diffusion or a validated local diffusivity method as appropriate.

## Production requirements

Provide a relaxed periodic membrane slab, validated bonded topology and charges, explicitly parameterized PFAS anions, compatible water/ion parameters, equilibrated starting structures and compute resources. GROMACS is not installed in the build environment. The supplied settings are preparation templates, not a validated force field or a completed simulation. Run and audit the protocol before reporting physical findings.

Code is MIT licensed; source data retains upstream attribution and terms. Large trajectories are excluded from Git and should be archived with hashes and persistent dataset identifiers.
