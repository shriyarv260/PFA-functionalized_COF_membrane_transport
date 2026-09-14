# Data contract

Every system must have a unique ID and a metadata JSON containing: species/CID, simulated protonation and net charge; COF ID and source checksum; functionalization sites and density; force-field versions and parameter hashes; water/ion models; box dimensions, membrane thickness, accessible area and pore measurement; solute concentration; temperature and ionic strength; seed; engine/version; time step; equilibration and production durations; restraint definitions; and raw-file SHA256 hashes.

Trajectory input: `positions_nm` is an unwrapped membrane-referenced molecular COM array (frames, molecules, 3), not wrapped atom coordinates. `dt_ps` is a positive constant output interval. Keep molecule ordering constant. MDAnalysis uses Angstrom coordinates by default, so divide by ten for nm. Validate the conversion against an engine distance output.

Umbrella input: `path,center_nm,k_kJ_mol_nm2`, with one finite numeric coordinate per line in each file. A GROMACS XVG contains time and coordinate columns and metadata; export only the intended coordinate column after equilibration/subsampling. Do not feed force values or time as positions. Histograms must include every retained sample. Specify bulk reference intervals separately when computing barriers from PMFs; the current CLI shifts the global minimum to zero.

Report per-replica estimates before aggregation. Use missing values for unavailable physical results. Never fill missing measured transport values with the synthetic demo. `results/public/` is descriptive database analysis; `results/synthetic/` is mathematical verification; `results/real/` is reserved for real simulation analyses.
