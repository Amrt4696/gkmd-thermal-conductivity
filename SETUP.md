# Setup

Building LAMMPS to run an MLIP like MACE is the annoying part - there's no single "pip install and go" here, since it depends on how you've exported your model and which pair style you're targeting.

## What this repo assumes

The LAMMPS input scripts here use:

```
pair_style mliap unified ${model} 0
pair_coeff * * <element list, in the order your model was trained on>
```

That means your LAMMPS build needs the ML-IAP package with the unified (Python-backed) interface enabled, and your MACE model needs to already be exported into the format `mliap unified` expects. The MACE -> LAMMPS export step itself lives in whatever repo you trained the potential in - out of scope here, this repo starts from an already-trained `.pt` model.

## CPU build

- cmake build with `-DPKG_ML-IAP=yes -DMLIAP_ENABLE_PYTHON=yes`, plus whichever other packages your workflow needs (KSPACE if you use long-range electrostatics elsewhere, etc.)
- needs a Python environment with `torch` and your model's dependencies importable from the same environment LAMMPS is linked against

## GPU build

- same ML-IAP setup, plus Kokkos/GPU package support depending on your LAMMPS version, and a CUDA-enabled torch matching your driver
- build with `-DBUILD_MPI=ON` even if you're only ever running single-GPU jobs - an MPI-disabled build is a common, silent cause of multi-GPU failures if you ever need to scale a system across cards later
- if you're running many independent replicas (different seeds, different P-T points) rather than one large system, it's usually simpler to run each as its own single-GPU process than to split one system's domain across multiple GPUs via MPI

## Environment notes

- `OMP_NUM_THREADS` matters more than you'd expect on a shared node - worth a short benchmark before committing to a value
- keep a separate build for anything you're actively rebuilding/testing vs. the binary your production jobs are pointed at, so a bad rebuild doesn't take down a run that's mid-flight

## Sanity check before trusting a build

Run a small system (a few hundred atoms) for a few thousand steps in NVE and watch `etotal` in the log. If it's drifting, something's off with the pair style or model export - worth catching here rather than 500 ps into a production run.
