"""Generate GROMACS window settings, never invent coordinates or parameters."""
import argparse
from pathlib import Path
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', required=True)
    p.add_argument('--range', nargs=2, type=float, default=[-2.,2.])
    p.add_argument('--spacing', type=float, default=.1)
    p.add_argument('--spring', type=float, default=1000.)
    p.add_argument('--membrane-pbcatom', type=int, required=True)
    p.add_argument('--pfas-pbcatom', type=int, required=True)
    a = p.parse_args()
    if a.spacing <= 0 or a.spring <= 0 or a.range[0] >= a.range[1] or min(a.membrane_pbcatom,a.pfas_pbcatom) < 1:
        p.error('Invalid range, spacing, spring or 1-based PBC reference atom')
    steps = round((a.range[1]-a.range[0])/a.spacing)
    if not np.isclose(steps*a.spacing,a.range[1]-a.range[0]): p.error('Spacing must divide range')
    base = (Path(__file__).resolve().parents[1]/'simulation/base.mdp').read_text()
    out = Path(a.out); out.mkdir(parents=True,exist_ok=True)
    for i,z in enumerate(np.linspace(*a.range,steps+1)):
        folder = out/f'window_{i:03d}'; folder.mkdir(exist_ok=True)
        pull = f'''
pull = yes
pull-ngroups = 2
pull-ncoords = 1
pull-group1-name = MEMBRANE
pull-group2-name = PFAS
pull-group1-pbcatom = {a.membrane_pbcatom}
pull-group2-pbcatom = {a.pfas_pbcatom}
pull-coord1-type = umbrella
pull-coord1-geometry = direction
pull-coord1-groups = 1 2
pull-coord1-dim = N N Y
pull-coord1-vec = 0 0 1
pull-coord1-start = no
pull-coord1-init = {z:.6f}
pull-coord1-rate = 0
pull-coord1-k = {a.spring}
pull-nstxout = 100
pull-nstfout = 100
'''
        (folder/'umbrella.mdp').write_text(base+pull)
    print(f'Created {steps+1} window templates. Supply physically positioned, equilibrated start.gro and matching topology/index files per window.')


if __name__ == '__main__': main()
