## bridgmanite_ambient_to_lowermantle_v1

MACE interatomic potential for bridgmanite (MgSiO3), trained on AIMD data spanning ambient conditions up to lower-mantle pressure-temperature conditions (135 GPa, 4000 K), across a broad set of configurations rather than a single P-T point.

Evaluated on an independent test set (structures held out of training):

| quantity | RMSE |
|---|---|
| energy  | 5.84 meV/atom |
| forces  | 7.71 meV/Å |
| stress  | 0.30 meV/Å^3 |

Exported for LAMMPS via the `mliap unified` interface - see `SETUP.md` for what your LAMMPS build needs to load it.

If you use this potential, a citation back to this repo (and the paper, once it's out) is appreciated.
