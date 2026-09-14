# Primary sources and roles

1. [CURATED-COFs](https://github.com/danieleongari/CURATED-COFs): framework registry, original structures, modifications and source papers. Ongari et al., ACS Central Science (2019), [doi:10.1021/acscentsci.9b00619](https://doi.org/10.1021/acscentsci.9b00619). Downloaded registry and license are pinned to a commit. Inspect geometry before simulation.
2. [Materials Cloud CURATED COFs](https://www.materialscloud.org/discover/curated-cofs): optimized structures and charge data referenced by CURATED. Selection and validation remain required; no optimized structure is claimed downloaded here.
3. [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest): molecular identities and properties for seven named parent acids. Raw responses, retrieval dates, hashes and CIDs are included. These are not force-field parameters.
4. [GROMACS umbrella-sampling tutorial](https://tutorials.gromacs.org/docs/umbrella-sampling.html): preparation and sampling methodology.
5. [GROMACS gmx wham documentation](https://manual.gromacs.org/current/onlinehelp/gmx-wham.html): PMF reconstruction, histograms, autocorrelation and bootstrap options. Save the exact engine manual/version with a production run.
6. [GROMACS MDP options](https://manual.gromacs.org/current/user-guide/mdp-options.html): parameter definitions and pull-coordinate conventions.

The public sources do not supply this project's historical MD trajectories or a validated PFAS–functionalized-COF transport benchmark. The database assessment therefore stops at identities, coverage and preparation requirements. Upstream materials carry their own licenses/terms; the repository MIT license applies to original code, not a blanket relicensing of external data.
