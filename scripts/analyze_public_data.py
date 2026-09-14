"""Descriptive statistics only: structural registries contain no transport labels."""
import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    folder = ROOT/'data/public'
    with (folder/'pfas_properties.csv').open() as f: pfas = list(csv.DictReader(f))
    with (folder/'cofs/cof-frameworks.csv').open() as f: cofs = list(csv.DictReader(f))
    summary = {'pfas_records':len(pfas), 'cof_records':len(cofs),
               'cof_columns':list(cofs[0]),
               'cof_nonempty_values':{k:sum(bool(r[k].strip()) for r in cofs) for k in cofs[0]},
               'transport_observations':0,
               'scope':'Downloaded registry metadata only; no historical project outputs supplied.'}
    out = ROOT/'results/public'; out.mkdir(parents=True, exist_ok=True)
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    lines = ['# Public data assessment', '',
             f'Analyzed {len(pfas)} PubChem parent-acid records and {len(cofs)} CURATED-COF registry rows.', '',
             '| Species | CID | Formula (parent acid) | Molecular weight (g/mol) |',
             '|---|---:|---|---:|']
    for r in sorted(pfas, key=lambda r:float(r['MolecularWeight'])):
        lines.append(f"| {r['species']} | {r['CID']} | {r['MolecularFormula']} | {r['MolecularWeight']} |")
    lines += ['', 'The selected panel spans carboxylic and sulfonic headgroups and multiple chain lengths. These descriptors support a controlled simulation design; they do not demonstrate rejection or adsorption trends.', '',
              '## Missing evidence', '',
              'No transport trajectories, force-field files, umbrella windows, pressure-driven crossing counts, or experimental rejection measurements were supplied. Consequently there are zero measured diffusion coefficients, PMFs, residence-time distributions, or water-permeance values in this dataset.', '',
              '## Structure triage', '',
              'Registry entries are candidate crystal structures, not prepared membrane systems. The upstream project notes possible geometry problems and points to optimized structures on Materials Cloud. Check atom overlaps, bonds, stacking, aqueous stability, accessible pore diameter, functionalization sites, charge, and finite-slab termination before selecting candidates.', '',
              '## Provenance', '',
              'See `data/public/manifest.json` for exact source URLs, retrieval times and SHA256 hashes; the COF revision is fixed in `data/sources.lock.json`. PubChem properties describe the retrieved neutral parent acids, not the proposed -1 simulation species.']
    (out/'assessment.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__': main()
