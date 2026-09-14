# Public data assessment

Analyzed 7 PubChem parent-acid records and 874 CURATED-COF registry rows.

| Species | CID | Formula (parent acid) | Molecular weight (g/mol) |
|---|---:|---|---:|
| PFBA | 9777 | C4HF7O2 | 214.04 |
| PFBS | 67815 | C4HF9O3S | 300.10 |
| PFHxA | 67542 | C6HF11O2 | 314.05 |
| PFHxS | 67734 | C6HF13O3S | 400.12 |
| PFOA | 9554 | C8HF15O2 | 414.07 |
| PFNA | 67821 | C9HF17O2 | 464.08 |
| PFOS | 74483 | C8HF17O3S | 500.13 |

The selected panel spans carboxylic and sulfonic headgroups and multiple chain lengths. These descriptors support a controlled simulation design; they do not demonstrate rejection or adsorption trends.

## Missing evidence

No transport trajectories, force-field files, umbrella windows, pressure-driven crossing counts, or experimental rejection measurements were supplied. Consequently there are zero measured diffusion coefficients, PMFs, residence-time distributions, or water-permeance values in this dataset.

## Structure triage

Registry entries are candidate crystal structures, not prepared membrane systems. The upstream project notes possible geometry problems and points to optimized structures on Materials Cloud. Check atom overlaps, bonds, stacking, aqueous stability, accessible pore diameter, functionalization sites, charge, and finite-slab termination before selecting candidates.

## Provenance

See `data/public/manifest.json` for exact source URLs, retrieval times and SHA256 hashes; the COF revision is fixed in `data/sources.lock.json`. PubChem properties describe the retrieved neutral parent acids, not the proposed -1 simulation species.
