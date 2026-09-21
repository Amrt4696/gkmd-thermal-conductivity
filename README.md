# GKMD Thermal Conductivity

Green-Kubo molecular dynamics workflow for computing lattice thermal conductivity of minerals, using LAMMPS driven by a machine-learned interatomic potential (MACE) instead of a classical force field.

I put this together while computing thermal conductivity of lower-mantle minerals (bridgmanite, MgO) for my PhD, and kept running into the same gap other people in this space seem to hit: there's plenty written about the Green-Kubo method in papers, and plenty about MLIPs, but not much showing the actual mechanics of wiring the two together in LAMMPS - equilibration, heat flux sampling, and the post-processing that turns a heat-flux time series into a kappa value with a real uncertainty estimate.

This repo is that workflow, generalized from what I actually run on our cluster. The worked example is bridgmanite (MgSiO3) at 120 GPa and 3500 K - roughly D'' layer conditions at the base of the mantle - using a MACE potential trained on our own AIMD data. Swap in your own structure, potential, and P-T point and the pipeline carries over.

## What's here

- `lammps/equilibration/` - three-stage cell relaxation (NPT, NPT, NVT) to get a well-equilibrated starting configuration at the target P-T
- `lammps/production/` - NVT production run with heat flux sampled every timestep, feeding the Green-Kubo integral
- `postprocess/kappa_from_heatflux.py` - FFT-based heat flux autocorrelation, Green-Kubo running integral, plateau averaging, anisotropy sanity check
- `docs/gotchas.md` - the debugging lessons that actually cost me time
- `SETUP.md` - getting LAMMPS built with an MLIP pair style in the first place, which is the part nobody documents well

## Method, briefly

Green-Kubo relates thermal conductivity to the time-autocorrelation of the heat flux:

```
kappa = V / (kB * T^2) * integral_0^inf <J(0) . J(t)> dt
```

In practice: equilibrate at the P-T point of interest, run NVT production while dumping the heat flux components every step, compute the autocorrelation function (FFT-based, since a direct sum over a million-step trajectory is too slow), and integrate it. The running integral - kappa as a function of the upper integration limit - climbs and then flattens into a plateau; you average over that plateau window to get kappa, and check that kxx/kyy/kzz agree with each other as a rough convergence check. If they don't, the run isn't long enough or the cell is too small.

For bridgmanite at 120 GPa / 3500 K, 10240 atoms, 1 ns production:

```
kappa = 6.52 +/- 0.03 W/m-K   (anisotropy ratio 1.19 across kxx/kyy/kzz)
```

## Quick start

1. Build LAMMPS with an MLIP pair style pointed at your MACE model - see `SETUP.md`.
2. Edit the P-T, cell, and model path at the top of `lammps/equilibration/input.lammps`, run it.
3. Point `lammps/production/input.lammps` at the resulting restart file, run it (expect >=1 ns; check the plateau before trusting the number).
4. `python postprocess/kappa_from_heatflux.py --file heatflux_*.dat --label mysystem --temperature 3500 --volume <cell volume in A^3>`

## Caveats

This is a research workflow, not a polished package - directory layout and defaults are stripped down for readability, not battle-tested against arbitrary systems or LAMMPS versions. If you hit a bug or run it on a system where something breaks, open an issue.
