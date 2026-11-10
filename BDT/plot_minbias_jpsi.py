#!/usr/bin/env python3
import argparse
import glob
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import uproot

from plot_minbias_no_bplus import BRANCHES, DATA, ENTRIES_PER_TASK, LOOPER, MINBIAS, MINBIAS_WEIGHT, OUT_DIR, STEP, \
    TOTAL_LOOPER_EVENTS, WORKERS, event_key, read_looper

CACHE = OUT_DIR / "minbias_jpsi.npz"
MASS_EDGES = np.linspace(2.0, 4.5, 126)
LXY_EDGES = np.array([0.0, 0.2, 1.0, 3.1, 7.0, 11.0, 16.0, 30.0, 70.0])
PLOT_BINS = (("0p0to0p2", r"$l_{xy}$ 0.0–0.2 cm", 0, 1), ("0p2to3p1", r"$l_{xy}$ 0.2–3.1 cm", 1, 3),
             ("3p1to11p0", r"$l_{xy}$ 3.1–11.0 cm", 3, 5), ("11p0to70p0", r"$l_{xy}$ 11.0–70.0 cm", 5, 8))
PEAK = (2.91, 3.27)
SIDEBANDS = ((2.65, 2.85), (3.31, 3.45))

KEYS = None
HAS_BPLUS = None
LIFETIME_WEIGHT = None


def fill(hist, mass, lxy, weight=1):
    ilxy = np.clip(np.digitize(lxy, LXY_EDGES) - 1, 0, len(LXY_EDGES) - 2)
    imass = np.clip(np.digitize(mass, MASS_EDGES) - 1, 0, len(MASS_EDGES) - 2)
    np.add.at(hist, (ilxy, imass), weight)


def scan(task):
    path, start, stop, is_mc = task
    shape = (len(LXY_EDGES) - 1, len(MASS_EDGES) - 1)
    total, no_bplus = np.zeros(shape, np.int64), np.zeros(shape, np.int64)
    reweighted, reweighted2 = np.zeros(shape), np.zeros(shape)
    with uproot.open(path) as f:
        for c in f["tuples"].iterate(BRANCHES, library="np", step_size=STEP, entry_start=start, entry_stop=stop):
            m, lxy = c["SV1_mass"], c["SV1_lxy"]
            hits = c["SV1_mu1_nhitsbeforesv"] + c["SV1_mu2_nhitsbeforesv"]
            hit_ok = np.where(lxy < 11.0, hits <= 0, np.where(lxy < 16.0, hits <= 1, hits <= 2))
            sel = hit_ok & (lxy >= 0.0) & (lxy <= 70.0) & (m >= MASS_EDGES[0]) & (m < MASS_EDGES[-1])
            if not is_mc:
                fill(total, m[sel], lxy[sel])
                continue
            key = event_key(c["lumi"], c["evtn"])
            i = np.clip(np.searchsorted(KEYS, key), 0, len(KEYS) - 1)
            sel &= KEYS[i] == key
            fill(total, m[sel], lxy[sel])
            keep = sel & ~HAS_BPLUS[i]
            fill(no_bplus, m[keep], lxy[keep])
            w = LIFETIME_WEIGHT[i[sel]].astype(np.float64)
            fill(reweighted, m[sel], lxy[sel], w)
            fill(reweighted2, m[sel], lxy[sel], w * w)
    return total, no_bplus, reweighted, reweighted2


def scan_all():
    global KEYS, HAS_BPLUS, LIFETIME_WEIGHT
    files = sorted(glob.glob(LOOPER))
    with ProcessPoolExecutor(8) as ex:
        looper = list(ex.map(read_looper, files))
    keys = np.concatenate([r[1] for r in looper])
    order = np.argsort(keys)
    KEYS, HAS_BPLUS = keys[order], np.concatenate([r[2] for r in looper])[order]
    LIFETIME_WEIGHT = np.concatenate([r[3] for r in looper])[order]
    counts = sum(r[0] for r in looper)
    print(f"{len(files)} looper files, {len(KEYS)} events, {counts / TOTAL_LOOPER_EVENTS:.1%} of MinBias", flush=True)
    with uproot.open(MINBIAS) as f:
        n = f["tuples"].num_entries
    with ProcessPoolExecutor(WORKERS) as ex:
        mc = list(ex.map(scan, [(MINBIAS, i, min(i + ENTRIES_PER_TASK, n), True)
                                for i in range(0, n, ENTRIES_PER_TASK)]))
        data = list(ex.map(scan, [(p, None, None, False) for p in sorted(glob.glob(DATA))], chunksize=4))
    return (sum(r[0] for r in mc), sum(r[1] for r in mc), sum(r[2] for r in mc), sum(r[3] for r in mc),
            sum(r[0] for r in data), counts)


def peak_yield(h, h2):
    lo, hi = MASS_EDGES[:-1], MASS_EDGES[1:]
    centers = 0.5 * (lo + hi)
    in_peak = (centers > PEAK[0]) & (centers < PEAK[1])
    in_side = np.zeros(len(centers), bool)
    for a, b in SIDEBANDS:
        in_side |= (centers > a) & (centers < b)
    scale = in_peak.sum() / in_side.sum()
    y = h[:, in_peak].sum(axis=1) - scale * h[:, in_side].sum(axis=1)
    v = h2[:, in_peak].sum(axis=1) + scale ** 2 * h2[:, in_side].sum(axis=1)
    return y, np.sqrt(v)


def plot_mass(mc_all, mc_veto, mc_rew, data, weight, show_reweighted=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep

    hep.style.use("CMS")
    centers = 0.5 * (MASS_EDGES[1:] + MASS_EDGES[:-1])
    for tag, label, i0, i1 in PLOT_BINS:
        d = data[i0:i1].sum(axis=0).astype(float)
        a = mc_all[i0:i1].sum(axis=0) * weight
        v = mc_veto[i0:i1].sum(axis=0) * weight
        r = mc_rew[i0:i1].sum(axis=0) * weight
        fig, (ax, rx) = plt.subplots(2, 1, figsize=(10, 9), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06})
        ax.stairs(a, MASS_EDGES, color="black", lw=1.5, ls="--", label=f"MinBias, all: {a.sum():,.0f}")
        ax.stairs(a - v, MASS_EDGES, color="#d62728", lw=2, label=f"MinBias, matched (gen $B^\\pm$ in event): {(a - v).sum():,.0f}")
        ax.stairs(v, MASS_EDGES, fill=True, alpha=0.35, color="tab:blue", label=f"MinBias, anti-matched (no gen $B^\\pm$): {v.sum():,.0f}")
        ax.stairs(v, MASS_EDGES, color="tab:blue")
        if show_reweighted:
            ax.stairs(r, MASS_EDGES, color="tab:green", lw=2, label=f"MinBias, $B^\\pm$ lifetime reweighted: {r.sum():,.0f}")
        nz = d > 0
        ax.errorbar(centers[nz], d[nz], yerr=np.sqrt(d[nz]), fmt="o", color="black", ms=3, label=f"Data: {d.sum():,.0f}")
        ax.set_yscale("log")
        ax.set_ylim(max(0.5, 0.3 * min(d[nz].min(), v[v > 0].min())), max(d.max(), a.max()) * 300)
        ax.set_ylabel("Events / 20 MeV")
        ax.legend(loc="upper right", fontsize=15, title=label, title_fontsize=15, frameon=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            rx.plot(centers[nz], (a / d)[nz], "o", ms=3, mfc="none", color="black", label="All")
            rx.plot(centers[nz], (v / d)[nz], "o", ms=3, color="tab:blue", label="Anti-matched")
            if show_reweighted:
                rx.plot(centers[nz], (r / d)[nz], "s", ms=3, color="tab:green", label="Reweighted")
        rx.axhline(1, color="gray", ls="--", lw=1)
        rx.set_yscale("log")
        rx.set_ylim(0.01, 100)
        rx.set_ylabel("MinBias / data", fontsize=16)
        rx.legend(loc="upper right", fontsize=12, ncol=3, frameon=False)
        rx.set_xlim(MASS_EDGES[0], MASS_EDGES[-1])
        rx.set_xlabel(r"SV1 $m_{\mu\mu}$ [GeV]")
        name = "minbias_jpsi_mass" if show_reweighted else "minbias_jpsi_noreweight_mass"
        fig.savefig(OUT_DIR / f"{name}_lxy_{tag}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


def plot_yield(mc_all, mc_veto, mc_rew, mc_rew2, data, weight):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep

    hep.style.use("CMS")
    yd, ed = peak_yield(data.astype(float), data.astype(float))
    ya, ea = peak_yield(mc_all * weight, mc_all * weight ** 2)
    yv, ev = peak_yield(mc_veto * weight, mc_veto * weight ** 2)
    yr, er = peak_yield(mc_rew * weight, mc_rew2 * weight ** 2)
    x = 0.5 * (LXY_EDGES[1:] + LXY_EDGES[:-1])
    xe = 0.5 * (LXY_EDGES[1:] - LXY_EDGES[:-1])
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(10, 9), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06})
    series = (("MinBias, all", ya, ea, "black", "o"), ("MinBias, anti-matched", yv, ev, "tab:blue", "o"),
              ("MinBias, $B^\\pm$ lifetime reweighted", yr, er, "tab:green", "s"))
    ax.errorbar(x, yd, xerr=xe, yerr=ed, fmt="D", color="gray", ms=6, label="Data")
    for label, y, e, color, marker in series:
        ax.errorbar(x, y, xerr=xe, yerr=e, fmt=marker, color=color, ms=5, label=label)
        with np.errstate(divide="ignore", invalid="ignore"):
            rx.errorbar(x, y / yd, xerr=xe, yerr=np.abs(y / yd) * np.sqrt((e / y) ** 2 + (ed / yd) ** 2), fmt=marker, color=color, ms=4)
    ax.set_yscale("log")
    ax.set_ylabel("Sideband-subtracted $J/\\psi$ events")
    ax.legend(loc="upper right", fontsize=15, frameon=False)
    rx.axhline(1, color="gray", ls="--", lw=1)
    rx.set_yscale("log")
    rx.set_ylim(0.01, 1000)
    rx.set_ylabel("MinBias / data", fontsize=16)
    rx.set_xscale("symlog", linthresh=1.0)
    rx.set_xlim(0, 70)
    rx.set_xticks([0, 1, 3.1, 7, 11, 16, 30, 70])
    rx.set_xticklabels(["0", "1", "3.1", "7", "11", "16", "30", "70"])
    rx.minorticks_off()
    rx.set_xlabel(r"SV1 $l_{xy}$ [cm]")
    fig.savefig(OUT_DIR / "minbias_jpsi_yield_vs_lxy.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("lxy bin, data, all, anti-matched, reweighted, ratios all/anti/rew")
    for i in range(len(x)):
        print(f"{LXY_EDGES[i]:5.1f}-{LXY_EDGES[i + 1]:5.1f} {yd[i]:11.0f}+-{ed[i]:7.0f} {ya[i]:11.0f}+-{ea[i]:7.0f} {yv[i]:11.0f}+-{ev[i]:7.0f} "
              f"{yr[i]:11.0f}+-{er[i]:7.0f}  {ya[i] / yd[i]:7.3f} {yv[i] / yd[i]:7.3f} {yr[i] / yd[i]:7.3f}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replot", action="store_true")
    args = parser.parse_args()
    OUT_DIR.mkdir(exist_ok=True)
    if not args.replot:
        mc_all, mc_veto, mc_rew, mc_rew2, data, counts = scan_all()
        np.savez(CACHE, mc_all=mc_all, mc_veto=mc_veto, mc_rew=mc_rew, mc_rew2=mc_rew2, data=data, counts=counts)
    with np.load(CACHE) as z:
        weight = MINBIAS_WEIGHT / (float(z["counts"]) / TOTAL_LOOPER_EVENTS)
        for show_reweighted in (True, False):
            plot_mass(z["mc_all"], z["mc_veto"], z["mc_rew"], z["data"], weight, show_reweighted)
        plot_yield(z["mc_all"], z["mc_veto"], z["mc_rew"], z["mc_rew2"], z["data"], weight)
    print(f"wrote {OUT_DIR}/minbias_jpsi_*.png")


if __name__ == "__main__":
    main()
