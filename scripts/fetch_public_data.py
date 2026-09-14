"""Fetch actual public records; preserve raw bytes, URLs, times and SHA256 hashes."""
import csv
import hashlib
import io
import json
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/public'
SPECIES = {'PFBA': 'perfluorobutanoic acid', 'PFHxA': 'perfluorohexanoic acid',
           'PFOA': 'perfluorooctanoic acid', 'PFNA': 'perfluorononanoic acid',
           'PFBS': 'perfluorobutanesulfonic acid', 'PFHxS': 'perfluorohexanesulfonic acid',
           'PFOS': 'perfluorooctanesulfonic acid'}


def fetch(url):
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent':'pfas-cof-research/0.1'})
            with urllib.request.urlopen(req, timeout=60) as r: return r.read()
        except Exception:
            if attempt == 3: raise
            time.sleep(2**attempt)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    def save(name, url):
        raw = fetch(url)
        (OUT/name).parent.mkdir(parents=True, exist_ok=True)
        (OUT/name).write_bytes(raw)
        manifest.append({'path':name, 'url':url, 'retrieved_utc':datetime.now(timezone.utc).isoformat(),
                         'sha256':hashlib.sha256(raw).hexdigest(), 'bytes':len(raw)})
        (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        return raw
    repo = 'danieleongari/CURATED-COFs'
    lock = ROOT/'data/sources.lock.json'
    if lock.exists(): sha = json.loads(lock.read_text())['curated_cofs_commit']
    else:
        sha = json.loads(fetch(f'https://api.github.com/repos/{repo}/commits/master'))['sha']
        lock.write_text(json.dumps({'curated_cofs_commit':sha}, indent=2)+'\n')
    base = f'https://raw.githubusercontent.com/{repo}/{sha}/'
    for file in ['cof-frameworks.csv','cof-papers.csv','LICENSE']:
        save('cofs/'+file, base+file)
    rows = []
    for label, name in SPECIES.items():
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.parse.quote(name)+'/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES,Charge/JSON'
        raw = save('pfas/'+label+'.json', url)
        records = json.loads(raw)['PropertyTable']['Properties']
        if len(records) != 1: raise ValueError(f'Ambiguous identity: {name}')
        record = records[0]
        rows.append({'species':label, **record, 'simulation_charge':-1,
                     'state_note':'Source is parent acid; simulation anion requires explicit parameterization.'})
        time.sleep(0.3)
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/'pfas_properties.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f'Saved {len(manifest)} source files and {len(rows)} PFAS records; COF commit {sha}')


if __name__ == '__main__': main()
