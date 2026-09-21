# Gotchas

Debugging lessons from actually running this, in no particular order.

**`timestep` after `read_restart`, not before.** If you set `timestep` before `read_restart`, the value gets silently overwritten by whatever was stored in the restart file - no warning, no error. You won't notice until you're computing kappa from a heat flux that was sampled at the wrong rate. Set it after.

**NVT, not NVE, for production.** Went back and forth on this. NVE conserves energy in principle, which sounds like the cleaner choice for a transport property, but in practice you don't get perfect energy conservation from an MLIP over a million steps, and the drift pollutes the heat flux more than a weakly-coupled thermostat does. Settled on NVT with loose damping for production and haven't revisited it.

**Weak thermostat damping during production, tight during settling.** Same reasoning - during production you want the thermostat doing as little as possible to the dynamics you're sampling, but during settling you want it doing real work to get the system to temperature quickly.

**Check energy conservation before trusting a new P-T point.** Before running full production, run the equilibrated structure in NVE for a few thousand steps and look at `etotal` in the log. If it's drifting, something's off with the potential or timestep, and it's worth catching before spending GPU-hours on a 1 ns run.

**Independent replicas, not domain decomposition, across GPUs.** If you're running several independent trajectories (different seeds, different P-T points) and each one fits comfortably on a single GPU, running each as its own single-GPU process is simpler and faster than splitting one system's MPI domain across multiple GPUs. Only reach for real multi-GPU decomposition when a single system genuinely doesn't fit in one card's memory.

**Shared LAMMPS binary + concurrent SLURM jobs = bus errors.** If multiple jobs load the same `liblammps.so` off a shared filesystem path at the same time, you can get intermittent bus errors that look like a hardware fault and aren't. Fix was giving each job its own copy of the binary in a per-job scratch directory instead of pointing everything at one shared path.

**Pick the plateau window by eye first.** `kappa_from_heatflux.py` takes a `--plateau START END` argument, but don't just guess it - plot the running integral first (the script saves a PNG) and find where it's actually flat. Too early and you're still climbing; too late and FFT noise in the tail starts dragging the average around. The anisotropy check (kxx vs kyy vs kzz) is a decent sanity flag if you're unsure - large disagreement between directions usually means the window is wrong or the run isn't long enough, not that the material is really that anisotropic.

**1 ns is a starting point, not a guarantee.** Whether the running integral has actually plateaued by 1 ns depends on the system - check the plot before trusting the number, and extend the run if it's still trending up.
