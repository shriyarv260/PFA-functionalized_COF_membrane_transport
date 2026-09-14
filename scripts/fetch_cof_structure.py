"""Download a selected registry CIF at the locked source revision."""
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from fetch_public_data import fetch, ROOT

p = argparse.ArgumentParser()
p.add_argument('cof_id', help='Exact CURATED-COFs ID from downloaded registry')
a = p.parse_args()
with (ROOT/'data/public/cofs/cof-frameworks.csv').open() as f:
    ids = {r['CURATED-COFs ID'] for r in csv.DictReader(f)}
if a.cof_id not in ids: p.error('ID is not present in the downloaded registry')
sha = json.loads((ROOT/'data/sources.lock.json').read_text())['curated_cofs_commit']
url = f'https://raw.githubusercontent.com/danieleongari/CURATED-COFs/{sha}/cifs/{a.cof_id}.cif'
raw = fetch(url)
out = ROOT/'data/structures'; out.mkdir(parents=True, exist_ok=True)
(out/f'{a.cof_id}.cif').write_bytes(raw)
(out/f'{a.cof_id}.provenance.json').write_text(json.dumps({
    'url':url, 'sha256':hashlib.sha256(raw).hexdigest(),
    'retrieved_utc':datetime.now(timezone.utc).isoformat(),
    'status':'Original registry CIF; unvalidated, not a prepared membrane.'}, indent=2)+'\n')
print(out/f'{a.cof_id}.cif')
