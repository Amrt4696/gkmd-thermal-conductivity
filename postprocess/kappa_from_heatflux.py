#!/usr/bin/env python3
"""
Green-Kubo thermal conductivity from a LAMMPS heat flux dump.

Takes the heat flux time series LAMMPS writes during an NVT production run
(via `fix ave/time ... file heatflux.dat`), computes the heat flux
autocorrelation function, integrates it (Green-Kubo), and estimates kappa
by averaging the running integral over a plateau window you pick by eye -
see docs/gotchas.md for how to pick it.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# physical constants
KB = 1.380649e-23          # J/K
EV_TO_J = 1.602176634e-19
PS_TO_S = 1.0e-12
A_TO_M = 1.0e-10


def load_heatflux(files, dt_sample_ps):
    frames = []
    for f in files:
        path = Path(f)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {f}")
        df = pd.read_csv(
            path, comment="#", sep=r"\s+",
            names=["step", "Jx", "Jy", "Jz"], engine="python",
        )
        frames.append(df)
        print(f"{f}: {len(df)} samples")

    data = pd.concat(frames, ignore_index=True)
    data["time_ps"] = np.arange(len(data)) * dt_sample_ps
    print(f"\nTotal samples : {len(data)}")
    print(f"Total time    : {data['time_ps'].iloc[-1]:.3f} ps")
    return data


def autocorr_fft(x, max_lag=None):
    """Heat flux autocorrelation via FFT - a direct sum over a
    million-step trajectory is too slow, this is the standard shortcut."""
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)

    n = len(x)
    f = np.fft.fft(x, n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[:n].real
    acf /= np.arange(n, 0, -1)  # unbiased normalization

    if max_lag is not None:
        acf = acf[:max_lag]
    return acf


def green_kubo_kappa(data, T, V_A3, dt_sample_ps, corr_ps):
    dt_sample_s = dt_sample_ps * 1e-12
    max_lag = int(corr_ps / dt_sample_ps)
    print(f"max_lag = {max_lag} points = {max_lag * dt_sample_ps:.3f} ps")

    Jx, Jy, Jz = (data[c].to_numpy() for c in ("Jx", "Jy", "Jz"))
    acf_x = autocorr_fft(Jx, max_lag)
    acf_y = autocorr_fft(Jy, max_lag)
    acf_z = autocorr_fft(Jz, max_lag)
    acf_avg = (acf_x + acf_y + acf_z) / 3.0
    tau_ps = np.arange(max_lag) * dt_sample_ps

    V_m3 = V_A3 * A_TO_M ** 3
    conv = EV_TO_J / (PS_TO_S * A_TO_M ** 2)   # eV/(ps*A^2) -> W/m^2
    prefactor = V_m3 * conv ** 2 / (KB * T ** 2)

    kappa_x = prefactor * np.cumsum(acf_x) * dt_sample_s
    kappa_y = prefactor * np.cumsum(acf_y) * dt_sample_s
    kappa_z = prefactor * np.cumsum(acf_z) * dt_sample_s
    kappa_avg = prefactor * np.cumsum(acf_avg) * dt_sample_s

    return tau_ps, kappa_x, kappa_y, kappa_z, kappa_avg


def plateau_average(tau_ps, kappa_x, kappa_y, kappa_z, kappa_avg,
                     plateau_start, plateau_end):
    mask = (tau_ps >= plateau_start) & (tau_ps <= plateau_end)

    k_x, k_y, k_z = (np.mean(k[mask]) for k in (kappa_x, kappa_y, kappa_z))
    k_mean = np.mean(kappa_avg[mask])
    k_std = np.std(kappa_avg[mask])
    anisotropy = max(k_x, k_y, k_z) / (min(k_x, k_y, k_z) + 1e-10)

    print(f"\nPlateau: {plateau_start:.1f} - {plateau_end:.1f} ps")
    print(f"kxx = {k_x:.3f} W/mK")
    print(f"kyy = {k_y:.3f} W/mK")
    print(f"kzz = {k_z:.3f} W/mK")
    print(f"kappa = {k_mean:.3f} +/- {k_std:.3f} W/mK")
    print(f"Anisotropy = {anisotropy:.2f}")

    if anisotropy > 2.0:
        print("WARNING: large anisotropy - likely not converged.")
    elif anisotropy > 1.5:
        print("Borderline anisotropy - treat result cautiously.")
    else:
        print("Anisotropy looks acceptable.")

    return k_mean, k_std, anisotropy


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--file", nargs="+", required=True, help="heat flux .dat file(s)")
    ap.add_argument("--label", required=True)
    ap.add_argument("--temperature", type=float, required=True, help="K")
    ap.add_argument("--volume", type=float, required=True, help="cell volume, A^3")
    ap.add_argument("--dt-sample-ps", type=float, default=0.001)
    ap.add_argument("--corr-ps", type=float, default=5.0, help="ACF window to compute")
    ap.add_argument("--plateau", type=float, nargs=2, default=[3.0, 4.0],
                     metavar=("START_PS", "END_PS"),
                     help="window to average the running integral over - look at "
                          "the plot first and pick where it's actually flat")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    data = load_heatflux(args.file, args.dt_sample_ps)
    tau_ps, kx, ky, kz, kavg = green_kubo_kappa(
        data, args.temperature, args.volume, args.dt_sample_ps, args.corr_ps
    )
    plateau_average(tau_ps, kx, ky, kz, kavg, *args.plateau)

    out = pd.DataFrame({"time_ps": tau_ps, "kappa_avg_WmK": kavg})
    out_name = f"TC_vs_time_avg_{args.label}.csv"
    out.to_csv(out_name, index=False)
    print(f"\nSaved: {out_name}")

    if not args.no_plot:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        axes[0].plot(tau_ps, kx, label="kxx")
        axes[0].plot(tau_ps, ky, label="kyy")
        axes[0].plot(tau_ps, kz, label="kzz")
        axes[0].set_xlabel("Upper integration limit (ps)")
        axes[0].set_ylabel("kappa (W/mK)")
        axes[0].legend()

        axes[1].plot(tau_ps, kavg, lw=2, color="red")
        axes[1].axvspan(*args.plateau, alpha=0.2, label="plateau window")
        axes[1].set_xlabel("Upper integration limit (ps)")
        axes[1].set_ylabel("kappa_avg (W/mK)")
        axes[1].legend()

        fig.tight_layout()
        fig.savefig(f"kappa_running_integral_{args.label}.png", dpi=150)
        print(f"Saved: kappa_running_integral_{args.label}.png")


if __name__ == "__main__":
    main()
