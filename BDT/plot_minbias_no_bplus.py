#!/usr/bin/env python3
import argparse
import glob
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import uproot

BASE = "/ceph/cms/store/group/Run3Scouting"
MINBIAS = f"{BASE}/tuples_minbias/tuples_MinBias_Fil-DoubleMuOS43_2024_2024.root"
DATA = f"{BASE}/tuples_data_nhfix/tuples_Data[CDEHI]_2024_*.root"
LOOPER = f"{BASE}/looperOutput_MinBias2024_genB/output_MinBias_Fil-DoubleMuOS43_2024_2024_*.root"
OUT_DIR = Path(__file__).resolve().parent / "minbias_vs_data"
CACHE = OUT_DIR / "minbias_no_bplus.npz"
RESONANCES = [(0.41, 0.50), (0.51, 0.59), (0.731, 0.83), (0.96, 1.08), (2.91, 3.27), (3.47, 3.89),
              (8.99, 9.91), (9.64, 10.56), (9.90, 10.78), (81.0, 101.0)]
LUMI_FB = 4.343
TOTAL_LOOPER_EVENTS = 235745961
MINBIAS_WEIGHT = 1.051e7 / 8.31e9 * 1000.0 * LUMI_FB
MASS_EDGES = np.unique(np.round(np.geomspace(0.2, 120.0, 121) / 0.002).astype(int)) * 0.002
LXY_EDGES = np.array([0.0, 0.2, 3.1, 11.0, 70.0])
LXY_BINS = (("0p0to0p2", r"$l_{xy}$ 0.0–0.2 cm"), ("0p2to3p1", r"$l_{xy}$ 0.2–3.1 cm"),
            ("3p1to11p0", r"$l_{xy}$ 3.1–11.0 cm"), ("11p0to70p0", r"$l_{xy}$ 11.0–70.0 cm"))
WORKERS = 24
STEP = 2_000_000
ENTRIES_PER_TASK = 4_000_000
BRANCHES = ["SV1_mass", "SV1_lxy", "SV1_mu1_nhitsbeforesv", "SV1_mu2_nhitsbeforesv", "lumi", "evtn"]

KEYS = None
HAS_BPLUS = None
LIFETIME_WEIGHT = None


def event_key(lumi, evtn):
    return (lumi.astype(np.int64) << 32) | evtn.astype(np.int64)


def read_looper(path):
    with uproot.open(path) as f:
        counts = float(f["counts"].values(flow=True).sum())
        a = f["tout"].arrays(["lumi", "evtn", "nGenBplus", "wGenBplus"], library="np")
    return counts, event_key(a["lumi"], a["evtn"]), a["nGenBplus"] > 0, a["wGenBplus"]


def fill(hist, mass, lxy, weight=1):
    ilxy = np.clip(np.digitize(lxy, LXY_EDGES[1:-1], right=True), 0, len(LXY_BINS) - 1)
    imass = np.clip(np.digitize(mass, MASS_EDGES) - 1, 0, len(MASS_EDGES) - 2)
    np.add.at(hist, (ilxy, imass), weight)


def scan(task):
    path, start, stop, is_mc = task
    shape = (len(LXY_BINS), len(MASS_EDGES) - 1)
    total, no_bplus = np.zeros(shape, np.int64), np.zeros(shape, np.int64)
    reweighted, reweighted2 = np.zeros(shape), np.zeros(shape)
    with uproot.open(path) as f:
        for c in f["tuples"].iterate(BRANCHES, library="np", step_size=STEP, entry_start=start, entry_stop=stop):
            m, lxy = c["SV1_mass"], c["SV1_lxy"]
            veto = np.zeros(len(m), bool)
            for lo, hi in RESONANCES:
                veto |= (m > lo) & (m < hi)
            hits = c["SV1_mu1_nhitsbeforesv"] + c["SV1_mu2_nhitsbeforesv"]
            hit_ok = np.where(lxy < 11.0, hits <= 0, np.where(lxy < 16.0, hits <= 1, hits <= 2))
            sel = ~veto & hit_ok & (lxy >= 0.0) & (lxy <= 70.0) & (m >= MASS_EDGES[0]) & (m < MASS_EDGES[-1])
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
    print(f"{len(files)} looper files, {len(KEYS)} events, {counts / TOTAL_LOOPER_EVENTS:.1%} of MinBias")
    with uproot.open(MINBIAS) as f:
        n = f["tuples"].num_entries
    with ProcessPoolExecutor(WORKERS) as ex:
        mc = list(ex.map(scan, [(MINBIAS, i, min(i + ENTRIES_PER_TASK, n), True)
                                for i in range(0, n, ENTRIES_PER_TASK)]))
        data = list(ex.map(scan, [(p, None, None, False) for p in sorted(glob.glob(DATA))], chunksize=4))
    return (sum(r[0] for r in mc), sum(r[1] for r in mc), sum(r[2] for r in mc), sum(r[3] for r in mc),
            sum(r[0] for r in data), counts)


def plot(mc_all, mc_veto, mc_rew, data, counts, show_reweighted=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep

    hep.style.use("CMS")
    frac = counts / TOTAL_LOOPER_EVENTS
    weight = MINBIAS_WEIGHT / frac
    centers = np.sqrt(MASS_EDGES[1:] * MASS_EDGES[:-1])
    ticks = [0.3, 0.5, 1, 2, 3, 5, 10, 20, 50, 100]
    for il, (tag, label) in enumerate(LXY_BINS):
        fig, (ax, rx) = plt.subplots(2, 1, figsize=(10, 9), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06})
        d, a, v = data[il].astype(float), mc_all[il] * weight, mc_veto[il] * weight
        nz = d > 0
        b = a - v
        r = mc_rew[il] * weight
        ax.stairs(a, MASS_EDGES, color="black", lw=1.5, ls="--", label=f"MinBias, all: {a.sum():,.0f}")
        ax.stairs(b, MASS_EDGES, color="#d62728", lw=2, label=f"MinBias, matched (gen $B^\\pm$ in event): {b.sum():,.0f}")
        ax.stairs(v, MASS_EDGES, fill=True, alpha=0.35, color="tab:blue", label=f"MinBias, anti-matched (no gen $B^\\pm$): {v.sum():,.0f}")
        ax.stairs(v, MASS_EDGES, color="tab:blue")
        if show_reweighted:
            ax.stairs(r, MASS_EDGES, color="tab:green", lw=2, label=f"MinBias, $B^\\pm$ lifetime reweighted: {r.sum():,.0f}")
        ax.errorbar(centers[nz], d[nz], yerr=np.sqrt(d[nz]), fmt="o", color="black", ms=3, label=f"Data: {d.sum():,.0f}")
        ax.set_yscale("log")
        ax.set_ylim(0.5, max(d.max(), a.max(), 1.0) * 30000)
        ax.set_ylabel("Events / bin")
        ax.legend(loc="upper right", fontsize=15, title=label if frac > 0.999 else f"{label}, {frac:.0%} of MinBias scaled up", title_fontsize=15,
                  frameon=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            rx.plot(centers[nz], (a / d)[nz], "o", ms=3, mfc="none", color="black", label="All")
            rx.plot(centers[nz], (v / d)[nz], "o", ms=3, color="tab:blue", label="Anti-matched")
            if show_reweighted:
                rx.plot(centers[nz], (r / d)[nz], "s", ms=3, color="tab:green", label="Reweighted")
        rx.axhline(1, color="gray", ls="--", lw=1)
        if il < 2:
            rx.set_ylim(0.0, 2.0)
        else:
            rx.set_yscale("log")
            rx.set_ylim(0.01, 100)
        rx.set_ylabel("MinBias / data", fontsize=16)
        rx.legend(loc="upper right", fontsize=12, ncol=3, frameon=False)
        rx.set_xscale("log")
        rx.set_xlim(MASS_EDGES[0], MASS_EDGES[-1])
        rx.set_xticks(ticks)
        rx.set_xticklabels([str(t) for t in ticks])
        rx.minorticks_off()
        rx.set_xlabel(r"SV1 $m_{\mu\mu}$ [GeV]")
        name = "minbias_no_bplus" if show_reweighted else "minbias_no_bplus_noreweight"
        fig.savefig(OUT_DIR / f"{name}_lxy_{tag}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    print(f"wrote {OUT_DIR}/{name}_lxy_*.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replot", action="store_true")
    args = parser.parse_args()
    OUT_DIR.mkdir(exist_ok=True)
    if not args.replot:
        mc_all, mc_veto, mc_rew, mc_rew2, data, counts = scan_all()
        np.savez(CACHE, mc_all=mc_all, mc_veto=mc_veto, mc_rew=mc_rew, mc_rew2=mc_rew2, data=data, counts=counts)
    with np.load(CACHE) as z:
        for show_reweighted in (True, False):
            plot(z["mc_all"], z["mc_veto"], z["mc_rew"], z["data"], float(z["counts"]), show_reweighted)


if __name__ == "__main__":
    main()
