# GROMACS execution contract

These templates cannot run without validated physical inputs. `base.mdp` assumes an already equilibrated NVT periodic system with compatible nonbonded settings; adjust cutoffs/dispersion corrections to the force field and document the change. It is not a minimization or initial-equilibration protocol.

```bash
python scripts/prepare_umbrellas.py --out runs/pfoa_parent/umbrellas --range -2 2 --spacing .1 --membrane-pbcatom 100 --pfas-pbcatom 1001
```

The atom numbers are examples only: use the actual 1-based reference atoms in the corresponding pull groups. Provide `start.gro`, `topol.top` and included parameter files, `index.ndx` containing MEMBRANE and a single PFAS molecule, and a compatible equilibrated checkpoint in each window. Audit initial signed COM distances and the box before running:

```bash
gmx grompp -f umbrella.mdp -c start.gro -t equil.cpt -p topol.top -n index.ndx -o umbrella.tpr
gmx mdrun -deffnm umbrella -px pullx.xvg -pf pullf.xvg
```

Run each window in its directory. Do not suppress grompp warnings with `-maxwarn`. The templates have not been compiled with GROMACS in this environment. Check sampling convergence independently of software validation; see `docs/protocol.md` for WHAM and production gates.
