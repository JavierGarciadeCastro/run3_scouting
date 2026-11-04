#!/usr/bin/env python3
import argparse
import glob
import json
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

_HERE = Path(__file__).resolve().parent
BASE = Path("/ceph/cms/store/group/Run3Scouting")
DATA_GLOB = str(BASE / "tuples_data" / "tuples_Data*_2024.root")
MINBIAS_FILE = str(BASE / "tuples_minbias" / "tuples_MinBias_Fil-DoubleMuOS43_2024_2024.root")
RESONANCE_WINDOWS = [(0.41, 0.50), (0.51, 0.59), (0.731, 0.83), (0.96, 1.08), (2.91, 3.27), (3.47, 3.89),
                     (8.99, 9.91), (9.64, 10.56), (9.90, 10.78), (81.0, 101.0)]
LXY_ORDER = ["0p0to0p2", "0p2to3p1", "3p1to11p0", "11p0to70p0"]
LUMI_FB = 109.95
UNBLIND_FRAC = 0.1
ERA_LUMI_FB = {"C": 7.24, "D": 7.96, "E": 11.32, "H": 5.44, "I": 11.47}
DATA_LUMI_FB = UNBLIND_FRAC * sum(ERA_LUMI_FB.values())
MINBIAS_NGEN = 409318867
MINBIAS_XSEC_PB = 1.051e7 * (MINBIAS_NGEN / 8.31e9)
HIST_EDGES = np.arange(0.0, 120.0 + 1e-9, 0.002)
EXTRA = ["SV1_mass", "SV1_lxy", "SV1_mu1_pt", "SV1_mu2_pt", "passL1", "SV1_mu1_nhitsbeforesv", "SV1_mu2_nhitsbeforesv"]
WP_MODE = "tight"
KEEP_MIN_TIGHT = 0.9
WORKERS = 10
THREADS = 4
ENTRIES_PER_TASK = 10_000_000
STEP = 500_000
PLOT_RANGE_REL = 0.5
PLOT_BINS = 40
FULL_BINS = 120


def _wp_dir(ev):
    return ev["dir"] / ("data_at_wp" if WP_MODE == "tight" else f"data_at_wp_{WP_MODE}")


def _wp_text(r):
    return f"FPR {r['fpr_target']:g}" if WP_MODE == "tight" else f"{r['sig_eff']:.0%} sig. eff."


def _flabel(f):
    return f"{f:g}".replace(".", "p")


def load_eval(eval_dir):
    eval_dir = Path(eval_dir).resolve()
    with open(eval_dir / "majority_vote_wp.json") as fh:
        wp = json.load(fh)
    thetas = sorted({tuple(v) for e in wp["mass_points"].values() for v in e.get("cond_theta", {}).values()})
    return {"dir": eval_dir, "wp": wp, "model_dir": wp["model_dir"], "features": wp["features"],
            "cond_vars": wp["cond_vars"], "conditional": wp["conditional"],
            "thetas": thetas if wp["conditional"] else [None],
            "require_l1": wp["require_l1"],
            "hit_preselection": wp.get("hit_preselection", False),
            "lxy_edges": [wp["lxy_bins"][LXY_ORDER[0]][0]] + [wp["lxy_bins"][lb][1] for lb in LXY_ORDER]}


def scorer_key(ev, theta):
    return "nominal" if theta is None else "theta_" + "_".join(_flabel(v) for v in theta)


def _cond_cols(ev):
    return [f"param_{v}" if not v.startswith("param_") else v for v in ev["cond_vars"]]


def _branches(evals, available, is_mc):
    want = set(EXTRA)
    for ev in evals:
        for f in ev["features"]:
            if f in _cond_cols(ev):
                continue
            if f.endswith("_dxy_lxy"):
                sv = f[:3]
                want |= {f"{sv}_lxy", f"{sv}_mass", f"{sv}_ptmm", f.replace("_dxy_lxy", "_dxy")}
            else:
                want.add(f)
    if is_mc:
        want.add("evtn")
    missing = sorted(b for b in want if b not in available)
    if missing:
        raise SystemExit(f"input is missing branches {missing[:10]}")
    return sorted(want)


def _add_dxy_lxy(df):
    for sv in ("SV1", "SV2"):
        if f"{sv}_ptmm" not in df:
            continue
        den = df[f"{sv}_lxy"] * df[f"{sv}_mass"] / df[f"{sv}_ptmm"]
        den = np.where(den > 1e-9, den, np.float32(1e-9))
        for mu in ("mu1", "mu2"):
            if f"{sv}_{mu}_dxy" in df:
                df[f"{sv}_{mu}_dxy_lxy"] = np.abs(df[f"{sv}_{mu}_dxy"]) / den


def _veto_mask(m):
    hit = np.zeros(len(m), bool)
    for lo, hi in RESONANCE_WINDOWS:
        hit |= (m > lo) & (m < hi)
    return hit


def score_chunk_task(task):
    path, sample, start, stop, eval_dirs, keep_min, threads, step = task
    import torch
    torch.set_num_threads(threads)
    sys.path.insert(0, str(_HERE))
    from train_discriminant import PNNClassifier
    evals = [load_eval(d) for d in eval_dirs]
    models = {}
    for i, ev in enumerate(evals):
        for lb in LXY_ORDER:
            models[(i, lb)] = PNNClassifier.load(str(Path(ev["model_dir"]) / f"pnn_lxy_{lb}.pt"))
    is_mc = sample == "minbias"
    out = {"hist": np.zeros((len(LXY_ORDER), len(HIST_EDGES) - 1), np.int64),
           "n_all_veto": 0, "n_test_veto": 0, "kept": {}}
    with uproot.open(path) as f:
        t = f["tuples"]
        br = _branches(evals, set(t.keys()), is_mc)
        for ch in t.iterate(br, library="np", step_size=step, entry_start=start, entry_stop=stop):
            evtn = ch.pop("evtn", None)
            df = pd.DataFrame({k: v.astype(np.float32) for k, v in ch.items()})
            m = df["SV1_mass"].to_numpy()
            sel = ~_veto_mask(m)
            if evals[0]["require_l1"]:
                sel &= df["passL1"].to_numpy() > 0.5
            if evals[0]["hit_preselection"]:
                nh = df["SV1_mu1_nhitsbeforesv"].to_numpy() + df["SV1_mu2_nhitsbeforesv"].to_numpy()
                lx = df["SV1_lxy"].to_numpy()
                sel &= np.where(lx < 11.0, nh <= 0, np.where(lx < 16.0, nh <= 1, nh <= 2))
            if is_mc:
                out["n_all_veto"] += int(sel.sum())
                sel &= (evtn % 10) >= 7
                out["n_test_veto"] += int(sel.sum())
            lxy = df["SV1_lxy"].to_numpy()
            edges = evals[0]["lxy_edges"]
            sel &= (lxy >= edges[0]) & (lxy <= edges[-1])
            df = df[sel].reset_index(drop=True)
            if not len(df):
                continue
            _add_dxy_lxy(df)
            b = np.digitize(df["SV1_lxy"].to_numpy(), edges[1:-1], right=True)
            m = df["SV1_mass"].to_numpy()
            for ib in range(len(LXY_ORDER)):
                out["hist"][ib] += np.histogram(m[b == ib], HIST_EDGES)[0]
            for i, ev in enumerate(evals):
                for theta in ev["thetas"]:
                    key = (i, scorer_key(ev, theta))
                    if theta is not None:
                        for c, v in zip(_cond_cols(ev), theta):
                            df[c] = np.float32(v)
                    X = df[ev["features"]].to_numpy(dtype=np.float32)
                    sc = np.full(len(df), -1.0, np.float32)
                    for ib, lb in enumerate(LXY_ORDER):
                        rows = np.flatnonzero(b == ib)
                        if len(rows):
                            sc[rows] = models[(i, lb)].predict_proba(X[rows])[:, 1]
                    k = sc > keep_min
                    piece = {"mass": m[k], "lxy": df["SV1_lxy"].to_numpy()[k], "bin": b[k].astype(np.int8),
                             "score": sc[k], "mu1_pt": df["SV1_mu1_pt"].to_numpy()[k],
                             "mu2_pt": df["SV1_mu2_pt"].to_numpy()[k]}
                    out["kept"].setdefault(key, []).append(piece)
    for key, pcs in out["kept"].items():
        out["kept"][key] = {c: np.concatenate([p[c] for p in pcs]) for c in pcs[0]}
    print(f"[score] {Path(path).name} [{start}, {stop}) done", flush=True)
    return sample, out


def run_scoring(evals, args, keep_min):
    files = [(p, "data") for p in sorted(glob.glob(args.data))]
    files.append((MINBIAS_FILE, "minbias"))
    tasks = []
    for path, sample in files:
        with uproot.open(path) as f:
            n = f["tuples"].num_entries
        for s in range(0, n, ENTRIES_PER_TASK):
            tasks.append((path, sample, s, min(n, s + ENTRIES_PER_TASK), [str(e["dir"]) for e in evals],
                          keep_min, THREADS, STEP))
    print(f"[score] {len(tasks)} tasks over {len(files)} files, {WORKERS} workers x {THREADS} threads")
    acc = {s: {"hist": 0, "n_all_veto": 0, "n_test_veto": 0, "kept": {}} for s in ("data", "minbias")}
    with mp.get_context("spawn").Pool(WORKERS) as pool:
        for sample, out in pool.imap_unordered(score_chunk_task, tasks):
            a = acc[sample]
            a["hist"] = a["hist"] + out["hist"]
            a["n_all_veto"] += out["n_all_veto"]
            a["n_test_veto"] += out["n_test_veto"]
            for key, d in out["kept"].items():
                a["kept"].setdefault(key, []).append(d)
    for i, ev in enumerate(evals):
        save = {}
        for sample, a in acc.items():
            if isinstance(a["hist"], int):
                continue
            save[f"{sample}_hist"] = a["hist"]
            save[f"{sample}_n_all_veto"] = a["n_all_veto"]
            save[f"{sample}_n_test_veto"] = a["n_test_veto"]
            for (j, sk), pcs in a["kept"].items():
                if j != i:
                    continue
                for c in pcs[0]:
                    save[f"{sample}_{sk}_{c}"] = np.concatenate([p[c] for p in pcs])
        save["keep_min_score"] = keep_min
        out_dir = _wp_dir(ev)
        out_dir.mkdir(exist_ok=True)
        np.savez(out_dir / "scores.npz", **save)
        print(f"[score] wrote {out_dir / 'scores.npz'}")


def _intervals_minus(base, cuts):
    segs = [tuple(base)]
    for lo, hi in cuts:
        nxt = []
        for a, b in segs:
            if hi <= a or lo >= b:
                nxt.append((a, b))
                continue
            if lo > a:
                nxt.append((a, lo))
            if hi < b:
                nxt.append((hi, b))
        segs = nxt
    return segs


def _count_in(m, segs):
    k = np.zeros(len(m), bool)
    for a, b in segs:
        k |= (m >= a) & (m <= b)
    return int(k.sum())


def _hist_count(h, segs):
    c = 0
    centers = 0.5 * (HIST_EDGES[1:] + HIST_EDGES[:-1])
    for a, b in segs:
        c += int(h[(centers >= a) & (centers <= b)].sum())
    return c


def _point_threshold(be):
    if WP_MODE == "loose":
        ts = sorted(set(be["threshold_loose_by_ctau"].values()))
    else:
        ts = sorted({be["threshold"]} if "threshold" in be else set(be["threshold_by_ctau"].values()))
    return float(np.median(ts)), len(ts)


def _point_scorer(ev, e):
    thetas = {tuple(v) for v in e.get("cond_theta", {}).values()}
    theta = next(iter(thetas)) if thetas else None
    return scorer_key(ev, theta), theta


def analyse(ev, args):
    z = np.load(_wp_dir(ev) / "scores.npz")
    wp = ev["wp"]
    veto = RESONANCE_WINDOWS
    has_mb = "minbias_hist" in z.files
    mb_w = 0.0
    if has_mb and int(z["minbias_n_test_veto"]) > 0:
        test_frac = int(z["minbias_n_test_veto"]) / int(z["minbias_n_all_veto"])
        mb_w = MINBIAS_XSEC_PB * 1e3 * LUMI_FB / (MINBIAS_NGEN * test_frac)
    w10 = mb_w * DATA_LUMI_FB / LUMI_FB
    keep_min = float(z["keep_min_score"])
    results = {}
    lines = [f"{'point':<16} {'lxy bin':>11} {'FPR':>7} | {'data before':>11} {'data after':>10} | "
             f"{'MB before':>11} {'MB after':>8} {'MB raw':>6}   (window, data lumi {DATA_LUMI_FB:.2f} fb^-1)"]
    for tag, e in sorted(wp["mass_points"].items(), key=lambda kv: kv[1].get("mA", kv[1].get("mzd"))):
        lo, hi = e["mass_window"]
        win = _intervals_minus((lo, hi), veto)
        sk, theta = _point_scorer(ev, e)
        res = {"mass_window": [lo, hi], "theta": list(theta) if theta else None, "bins": {}}
        for ib, lb in enumerate(LXY_ORDER):
            be = e["bins"].get(lb)
            if be is None:
                continue
            thr, n_thr = _point_threshold(be)
            if n_thr > 1:
                print(f"[warn] {tag} {lb}: {n_thr} different thresholds across ctau, using the median")
            if thr <= keep_min:
                raise SystemExit(f"{tag} {lb}: threshold {thr} is below keep_min_score {keep_min}; rerun with --rescore or lower KEEP_MIN_TIGHT")
            r = {"fpr_target": be["fpr_target"], "threshold": thr, "wp": WP_MODE,
                 "sig_eff": wp.get("loose_sig_eff") if WP_MODE == "loose" else None}
            for sample in ("data",) + (("minbias",) if has_mb else ()):
                p = (z[f"{sample}_{sk}_bin"] == ib) & (z[f"{sample}_{sk}_score"] > thr)
                mp_ = z[f"{sample}_{sk}_mass"][p]
                h = z[f"{sample}_hist"][ib]
                n_bin = int(h.sum())
                d = {"n_window_before_wp": _hist_count(h, win), "n_window_after_wp": _count_in(mp_, win),
                     "n_bin_before_wp": n_bin, "n_bin_after_wp": int(len(mp_)),
                     "pass_frac_bin": len(mp_) / n_bin if n_bin else None}
                if sample == "minbias":
                    d["weight_per_event_data_lumi"] = w10
                    d["n_window_before_wp_data_lumi"] = d["n_window_before_wp"] * w10
                    d["n_window_after_wp_data_lumi"] = d["n_window_after_wp"] * w10
                r[sample] = d
            res["bins"][lb] = r
            dd, mm = r["data"], r.get("minbias")
            wp_col = f"{be['fpr_target']:7.0e}" if WP_MODE == "tight" else f"{thr:7.4f}"
            lines.append(f"{tag:<16} {lb:>11} {wp_col} | {dd['n_window_before_wp']:11,d} "
                         f"{dd['n_window_after_wp']:10d} | "
                         + (f"{mm['n_window_before_wp_data_lumi']:11,.0f} {mm['n_window_after_wp_data_lumi']:8.0f} "
                            f"{mm['n_window_after_wp']:6d}" if mm else ""))
        results[tag] = res
    out = {"description": "Data (10% of 2024 C,D,E,H,I) and MinBias test set in each signal mass window and lxy "
                          f"bin, before and after the majority-vote WP of {ev['dir']}. MinBias *_data_lumi = raw "
                          f"test-set counts x weight_per_event_data_lumi (eval normalisation scaled to the data "
                          f"lumi, {DATA_LUMI_FB:.3f} fb^-1 = 10% of eras C,D,E,H,I, i.e. comparable to the data). "
                          "pass_frac_bin = after/before over the whole lxy "
                          "bin (all masses).",
           "minbias_weight_per_test_event_full_lumi": mb_w, "mass_points": results}
    path = _wp_dir(ev) / "data_at_wp.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\n==== {ev['dir']}")
    print("\n".join(lines))
    print(f"wrote {path}")
    plot_dir = _wp_dir(ev) / "plots"
    plot_mass_spectra(ev, z, out, w10, plot_dir)
    plot_mass_before_after(ev, z, out, w10, plot_dir)
    plot_full_spectrum(ev, z, out, w10, plot_dir)
    if has_mb:
        plot_pass_fraction(ev, out, plot_dir)
        plot_window_yields(ev, out, plot_dir)
    print(f"plots in {plot_dir}")
    if args.write_scored:
        write_scored(ev, z, args)


def _plot_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _disc_label(ev):
    return "PNN" if ev["conditional"] else "NN"


def _mass_title(e):
    if "mzd" in e:
        return rf"$m_{{Z_d}} = {e['mzd']:g}$ GeV"
    return rf"$m_{{\pi_3}} = {e['mpi']:g}$ GeV, $m_{{A'}} = {e['mA']:g}$ GeV"


def plot_mass_spectra(ev, z, out, w10, plot_dir):
    plt = _plot_style()
    d = plot_dir / "mass_after_wp"
    d.mkdir(parents=True, exist_ok=True)
    wp = ev["wp"]
    has_mb = "minbias_hist" in z.files
    for tag, res in out["mass_points"].items():
        e = wp["mass_points"][tag]
        x0, x1 = res["mass_window"]
        edges = np.linspace(x0, x1, PLOT_BINS + 1)
        c = 0.5 * (edges[1:] + edges[:-1])
        sk, _ = _point_scorer(ev, e)
        fig, axes = plt.subplots(2, 2, figsize=(13, 9.5), sharex=True)
        for ib, (lb, ax) in enumerate(zip(LXY_ORDER, axes.flat)):
            r = res["bins"].get(lb)
            if r is None:
                ax.set_visible(False)
                continue
            thr = r["threshold"]
            ymax = 1.0
            if has_mb:
                p = (z[f"minbias_{sk}_bin"] == ib) & (z[f"minbias_{sk}_score"] > thr)
                hm, _ = np.histogram(z[f"minbias_{sk}_mass"][p], edges)
                ax.stairs(hm * w10, edges, fill=True, alpha=0.35, color="tab:blue",
                          label=f"MinBias (1 raw evt = {w10:.0f})")
                ax.stairs(hm * w10, edges, color="tab:blue")
                ymax = max(ymax, float((hm * w10).max()))
            p = (z[f"data_{sk}_bin"] == ib) & (z[f"data_{sk}_score"] > thr)
            hd, _ = np.histogram(z[f"data_{sk}_mass"][p], edges)
            nz = hd > 0
            ax.errorbar(c[nz], hd[nz], yerr=np.sqrt(hd[nz]), fmt="o", color="black", ms=4, label=f"Data (10% of C,D,E,H,I, {DATA_LUMI_FB:.2f} fb$^{{-1}}$)")
            ymax = max(ymax, float(hd.max()) if len(hd) else 0.0)
            ax.set_ylim(0, ymax * 1.9 + 1)
            dd = r["data"]
            txt = (rf"$l_{{xy}}$ {lb.replace('p', '.').replace('to', '–')} cm, {_wp_text(r)}" + "\n"
                   + f"In window, before -> after WP:\n"
                   + f"  Data:    {dd['n_window_before_wp']:,} -> {dd['n_window_after_wp']}")
            if has_mb:
                mm = r["minbias"]
                txt += (f"\n  MinBias: {mm['n_window_before_wp_data_lumi']:,.0f} -> "
                        f"{mm['n_window_after_wp_data_lumi']:.0f} ({mm['n_window_after_wp']} raw)")
            ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left", fontsize=9, family="monospace")
            ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.74), fontsize=8, framealpha=0.9)
            ax.set_ylabel(f"Events / {(x1 - x0) / PLOT_BINS * 1000:.0f} MeV", fontsize=10)
            ax.tick_params(direction="in", top=True, right=True)
        for ax in axes[-1]:
            ax.set_xlabel(r"SV1 $m_{\mu\mu}$ [GeV]")
        fig.suptitle(f"{_disc_label(ev)} majority-vote WP, {_mass_title(e)}: events passing the WP", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(d / f"mass_after_wp_{tag}.png", dpi=110)
        plt.close(fig)


def plot_mass_before_after(ev, z, out, w10, plot_dir):
    plt = _plot_style()
    d = plot_dir / "mass_before_after"
    d.mkdir(parents=True, exist_ok=True)
    wp = ev["wp"]
    has_mb = "minbias_hist" in z.files
    hstep = HIST_EDGES[1] - HIST_EDGES[0]
    hc = 0.5 * (HIST_EDGES[1:] + HIST_EDGES[:-1])
    for tag, res in out["mass_points"].items():
        e = wp["mass_points"][tag]
        mref = e.get("mA", e.get("mzd"))
        lo, hi = res["mass_window"]
        x0 = max(0.0, lo - PLOT_RANGE_REL * mref)
        x1 = hi + PLOT_RANGE_REL * mref
        bw = max(1, round((x1 - x0) / PLOT_BINS / hstep)) * hstep
        x0 = np.floor(x0 / hstep) * hstep
        edges = np.arange(x0, x1 + bw, bw)
        c = 0.5 * (edges[1:] + edges[:-1])
        sk, _ = _point_scorer(ev, e)
        for ib, lb in enumerate(LXY_ORDER):
            r = res["bins"].get(lb)
            if r is None:
                continue
            fig, axes = plt.subplots(2, 2, figsize=(14, 6.5), sharex=True,
                                     gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06, "wspace": 0.18})
            for col, stage in enumerate(("before", "after")):
                ax, rx = axes[0, col], axes[1, col]
                hists = {}
                for sample in ("data",) + (("minbias",) if has_mb else ()):
                    if stage == "before":
                        h, _ = np.histogram(hc, edges, weights=z[f"{sample}_hist"][ib])
                    else:
                        p = (z[f"{sample}_{sk}_bin"] == ib) & (z[f"{sample}_{sk}_score"] > r["threshold"])
                        h, _ = np.histogram(z[f"{sample}_{sk}_mass"][p], edges)
                    hists[sample] = h.astype(float) * (w10 if sample == "minbias" else 1.0)
                hd = hists["data"]
                nz = hd > 0
                if has_mb:
                    hm = hists["minbias"]
                    ax.stairs(hm, edges, fill=True, alpha=0.35, color="tab:blue",
                              label=f"MinBias: {hm.sum():,.0f} ({hm.sum() / w10:,.0f} raw)")
                    ax.stairs(hm, edges, color="tab:blue")
                    with np.errstate(divide="ignore", invalid="ignore"):
                        rr = hm / hd
                    rx.plot(c[nz], rr[nz], "o", ms=3, color="tab:blue")
                ax.errorbar(c[nz], hd[nz], yerr=np.sqrt(hd[nz]), fmt="o", color="black", ms=3,
                            label=f"Data: {hd.sum():,.0f}")
                ax.set_yscale("log")
                top = max(hd.max(), hists["minbias"].max() if has_mb else 0.0, 1.0)
                ax.set_ylim(0.5, top * 30)
                for x_ in (ax, rx):
                    for v in (lo, hi):
                        x_.axvline(v, color="#d62728", lw=1.0, ls="--")
                rx.axhline(1, color="gray", ls="--", lw=1)
                rx.set_yscale("log")
                rx.set_ylim(0.05, 100)
                rx.set_ylabel("MinBias / data", fontsize=8)
                rx.tick_params(labelsize=8, direction="in", which="both", top=True)
                ax.tick_params(labelsize=8, direction="in", which="both", top=True, right=True)
                ax.set_ylabel(f"Events / {bw * 1000:.0f} MeV", fontsize=9)
                head = "Before WP" if stage == "before" else f"After WP ({_wp_text(r)})"
                ax.legend(loc="upper right", fontsize=9, title=head, title_fontsize=10)
            for rx in axes[-1]:
                rx.set_xlabel(r"SV1 $m_{\mu\mu}$ [GeV]")
            fig.suptitle(f"{_mass_title(e)}, " + rf"$l_{{xy}}$ {lb.replace('p', '.').replace('to', '–')} cm",
                         fontsize=13)
            (d / tag).mkdir(exist_ok=True)
            fig.savefig(d / tag / f"mass_before_after_{tag}_lxy_{lb}.png", dpi=100, bbox_inches="tight")
            plt.close(fig)


def plot_full_spectrum(ev, z, out, w10, plot_dir):
    plt = _plot_style()
    d = plot_dir / "full_spectrum"
    d.mkdir(parents=True, exist_ok=True)
    wp = ev["wp"]
    samples = ("data", "minbias") if "minbias_hist" in z.files else ("data",)
    hstep = HIST_EDGES[1] - HIST_EDGES[0]
    hc = 0.5 * (HIST_EDGES[1:] + HIST_EDGES[:-1])
    edges = np.unique(np.round(np.geomspace(0.2, 120.0, FULL_BINS + 1) / hstep).astype(int)) * hstep
    c = np.sqrt(edges[1:] * edges[:-1])
    before = {s: z[f"{s}_hist"] for s in samples}
    before = {s: [np.histogram(hc, edges, weights=h[ib])[0].astype(float) for ib in range(len(LXY_ORDER))]
              for s, h in before.items()}
    kept = {}
    for tag, res in out["mass_points"].items():
        e = wp["mass_points"][tag]
        lo, hi = res["mass_window"]
        sk, _ = _point_scorer(ev, e)
        for s in samples:
            if (s, sk) not in kept:
                kept[(s, sk)] = {c_: z[f"{s}_{sk}_{c_}"] for c_ in ("bin", "score", "mass")}
        fig, axes = plt.subplots(2 * len(LXY_ORDER), 2, figsize=(14, 19), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1] * len(LXY_ORDER), "hspace": 0.1, "wspace": 0.16})
        for ib, lb in enumerate(LXY_ORDER):
            r = res["bins"].get(lb)
            if r is None:
                for col in range(2):
                    axes[2 * ib, col].set_visible(False)
                    axes[2 * ib + 1, col].set_visible(False)
                continue
            for col, stage in enumerate(("before", "after")):
                ax, rx = axes[2 * ib, col], axes[2 * ib + 1, col]
                hists = {}
                for s in samples:
                    if stage == "before":
                        h = before[s][ib]
                    else:
                        k = kept[(s, sk)]
                        p = (k["bin"] == ib) & (k["score"] > r["threshold"])
                        h = np.histogram(k["mass"][p], edges)[0].astype(float)
                    hists[s] = h * (w10 if s == "minbias" else 1.0)
                hd = hists["data"]
                nz = hd > 0
                if "minbias" in hists:
                    hm = hists["minbias"]
                    ax.stairs(hm, edges, fill=True, alpha=0.35, color="tab:blue",
                              label=f"MinBias: {hm.sum():,.0f} ({hm.sum() / w10:,.0f} raw)")
                    ax.stairs(hm, edges, color="tab:blue")
                    with np.errstate(divide="ignore", invalid="ignore"):
                        rr = hm / hd
                    rx.plot(c[nz], rr[nz], "o", ms=3, color="tab:blue")
                ax.errorbar(c[nz], hd[nz], yerr=np.sqrt(hd[nz]), fmt="o", color="black", ms=3,
                            label=f"Data: {hd.sum():,.0f}")
                ax.set_yscale("log")
                ax.set_ylim(0.5, max(max(h.max() for h in hists.values()), 1.0) * 3000)
                for x_ in (ax, rx):
                    for v in (lo, hi):
                        x_.axvline(v, color="#d62728", lw=1.0, ls="--")
                rx.axhline(1, color="gray", ls="--", lw=1)
                rx.set_yscale("log")
                rx.set_ylim(0.05, 100)
                rx.set_ylabel("MinBias / data", fontsize=8)
                rx.tick_params(labelsize=8, direction="in", which="both", top=True)
                ax.tick_params(labelsize=8, direction="in", which="both", top=True, right=True)
                ax.set_ylabel("Events / bin", fontsize=9)
                rng = rf"$l_{{xy}}$ {lb.replace('p', '.').replace('to', '–')} cm"
                head = f"{rng}, before WP" if stage == "before" else f"{rng}, after WP ({_wp_text(r)})"
                ax.legend(loc="upper right", fontsize=9, title=head, title_fontsize=10)
        for ax in axes[-1]:
            ax.set_xscale("log")
            ax.set_xlim(edges[0], edges[-1])
            ax.set_xticks([0.3, 0.5, 1, 2, 3, 5, 10, 20, 50, 100])
            ax.set_xticklabels(["0.3", "0.5", "1", "2", "3", "5", "10", "20", "50", "100"])
            ax.minorticks_off()
            ax.set_xlabel(r"SV1 $m_{\mu\mu}$ [GeV]")
        fig.suptitle(f"{_disc_label(ev)}, {_mass_title(e)}: full mass spectrum, data lumi {DATA_LUMI_FB:.2f} fb$^{{-1}}$",
                     fontsize=13, y=0.905)
        fig.savefig(d / f"full_spectrum_{tag}.png", dpi=100, bbox_inches="tight")
        plt.close(fig)


def plot_window_yields(ev, out, plot_dir):
    plt = _plot_style()
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullFormatter
    pts = ev["wp"]["mass_points"]
    tags = list(out["mass_points"])
    is_dqcd = "mA" in pts[tags[0]]

    def _x(t):
        return pts[t]["mA"] if is_dqcd else pts[t]["mzd"]

    def _grp(t):
        return int(round(pts[t]["mpi"] / pts[t]["mA"])) if is_dqcd else None

    groups = sorted({_grp(t) for t in tags})
    styles = [("-", "o"), ("--", "s")]
    series = [("MinBias before WP", "tab:blue", lambda r: r["minbias"]["n_window_before_wp_data_lumi"]),
              ("Data before WP", "black", lambda r: r["data"]["n_window_before_wp"]),
              ("Data after WP", "#d62728", lambda r: r["data"]["n_window_after_wp"])]
    floor = 0.5
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    for lb, ax in zip(LXY_ORDER, axes.flat):
        for i, g in enumerate(groups):
            ls, mk = styles[i % 2]
            sel = sorted((t for t in tags if _grp(t) == g and lb in out["mass_points"][t]["bins"]), key=_x)
            xs = np.array([_x(t) for t in sel])
            for _, color, get in series:
                ys = np.array([get(out["mass_points"][t]["bins"][lb]) for t in sel], dtype=float)
                zero = ys <= 0
                ax.plot(xs[~zero], ys[~zero], color=color, ls=ls, marker=mk, ms=5, lw=1.4)
                ax.plot(xs[zero], np.full(zero.sum(), floor), color=color, ls="none", marker="v", ms=7, mfc="none")
        ax.set_title(rf"$l_{{xy}}$ {lb.replace('p', '.').replace('to', '–')} cm", loc="left", fontsize=13)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(floor / 2, None)
        ax.set_ylabel(f"Events in mass window ({DATA_LUMI_FB:.2f} fb$^{{-1}}$)")
        ax.grid(alpha=0.3, which="both")
        ax.tick_params(direction="in", top=True, right=True, which="both")
    handles = [Line2D([], [], color=color, lw=2, label=lab) for lab, color, _ in series]
    if is_dqcd:
        handles += [Line2D([], [], color="gray", ls=styles[i % 2][0], marker=styles[i % 2][1],
                           label=rf"$m_{{A'}}/m_{{\pi_3}} = 1/{g}$") for i, g in enumerate(groups)]
    handles.append(Line2D([], [], color="gray", ls="none", marker="v", mfc="none", ms=7, label="0 events"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.96), fontsize=11, ncol=len(handles),
               frameon=False)
    for ax in axes[-1]:
        if is_dqcd:
            ax.set_xticks([0.3, 0.5, 1, 2, 3])
            ax.set_xticklabels(["0.3", "0.5", "1", "2", "3"])
            ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(r"$m_{A'}$ [GeV]" if is_dqcd else r"$m_{Z_d}$ [GeV]")
    head = "majority-vote WP" if WP_MODE == "tight" else f"loose WP ({ev['wp']['loose_sig_eff']:.0%} sig. eff.)"
    fig.suptitle(f"{_disc_label(ev)}, {head}: events in each signal mass window, resonance windows removed", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(plot_dir / "window_yields.png", dpi=130)
    plt.close(fig)


def plot_pass_fraction(ev, out, plot_dir):
    plt = _plot_style()
    wp = ev["wp"]
    tags = sorted(out["mass_points"], key=lambda t: (wp["mass_points"][t].get("mA", wp["mass_points"][t].get("mzd")),
                                                     wp["mass_points"][t].get("mpi", 0)))
    x = np.arange(len(tags))
    fig, axes = plt.subplots(4, 1, figsize=(12, 13), sharex=True)
    for lb, ax in zip(LXY_ORDER, axes):
        vals = {"data": [], "minbias": [], "fpr": []}
        for t in tags:
            r = out["mass_points"][t]["bins"].get(lb)
            vals["fpr"].append(r["fpr_target"] if r and WP_MODE == "tight" else np.nan)
            for s in ("data", "minbias"):
                v = r[s]["pass_frac_bin"] if r and r.get(s) else None
                vals[s].append(v if v else np.nan)
        ax.plot(x, vals["fpr"], "_", color="0.4", ms=22, mew=2, label="FPR target")
        ax.plot(x, vals["minbias"], "s", color="tab:blue", ms=7, label="MinBias test set")
        ax.plot(x, vals["data"], "o", color="black", ms=6, label="Data (10%)")
        ax.set_yscale("log")
        ax.set_ylabel("Fraction passing WP")
        ax.grid(alpha=0.3, which="both", axis="y")
        ax.text(0.01, 0.95, rf"$l_{{xy}}$ {lb.replace('p', '.').replace('to', '–')} cm", transform=ax.transAxes,
                va="top", fontsize=12)
        ax.tick_params(direction="in", top=True, right=True, which="both")
    axes[0].legend(loc="upper right", fontsize=10, ncol=3)
    axes[0].set_title(f"{_disc_label(ev)}: fraction of events in the lxy bin passing the WP (all masses)",
                      fontsize=12, loc="left")
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(tags, rotation=60, ha="right", fontsize=9)
    fig.tight_layout()
    fig.savefig(plot_dir / "pass_fraction_data_vs_minbias.png", dpi=110)
    plt.close(fig)


def write_scored(ev, z, args):
    out_dir = _wp_dir(ev) / "scored"
    out_dir.mkdir(exist_ok=True)
    for tag, e in ev["wp"]["mass_points"].items():
        sk, _ = _point_scorer(ev, e)
        b = z[f"data_{sk}_bin"]
        sc = z[f"data_{sk}_score"]
        m = z[f"data_{sk}_mass"]
        thr = np.full(len(b), np.inf, np.float32)
        for ib, lb in enumerate(LXY_ORDER):
            be = e["bins"].get(lb)
            if be is not None:
                thr[b == ib] = _point_threshold(be)[0]
        keep = sc > thr
        cols = {"SV1_mass": m[keep], "SV1_lxy": z[f"data_{sk}_lxy"][keep],
                "SV1_mu1_pt": z[f"data_{sk}_mu1_pt"][keep], "SV1_mu2_pt": z[f"data_{sk}_mu2_pt"][keep],
                "bdt_score": sc[keep], "bdt_threshold": thr[keep],
                "lxy_bin_idx": b[keep].astype(np.float32), "pass_bdt": np.ones(int(keep.sum()), np.float32)}
        path = out_dir / f"scored_Data2024_{tag}.root"
        with uproot.recreate(str(path)) as f:
            t = f.mktree("tuples", {c: v.dtype for c, v in cols.items()})
            t.extend(cols)
        print(f"[scored] {path.name}: {int(keep.sum())} WP-passing data events")


def main():
    p = argparse.ArgumentParser(description="Score 2024 data (and the MinBias test set) with the NN/PNN models of "
                                            "one or more evaluate_workingpoint.py eval dirs and count events "
                                            "surviving each majority-vote WP.")
    p.add_argument("eval_dirs", nargs="+", help="eval_testonly* dirs holding majority_vote_wp.json")
    p.add_argument("--data", default=DATA_GLOB)
    p.add_argument("--wp", choices=["tight", "loose"], default="tight",
                   help="tight: majority-vote WP (threshold_by_ctau). loose: threshold_loose_by_ctau (loose_sig_eff "
                        "signal efficiency, from evaluate_workingpoint.py). Outputs go to data_at_wp/ or "
                        "data_at_wp_loose/.")
    p.add_argument("--rescore", action="store_true", help="Rescore even if scores.npz exists.")
    p.add_argument("--write-scored", action="store_true",
                   help="Write per-mass-point ROOT tuples of WP-passing data for make_fit_inputs.py --from-scored.")
    args = p.parse_args()
    global WP_MODE
    WP_MODE = args.wp

    evals = [load_eval(d) for d in args.eval_dirs]
    for ev in evals[1:]:
        for k in ("require_l1", "hit_preselection", "lxy_edges"):
            if ev[k] != evals[0][k]:
                raise SystemExit(f"eval dirs disagree on {k}; run them separately")
    if WP_MODE == "tight":
        keep_min = KEEP_MIN_TIGHT
    else:
        thrs = []
        for ev in evals:
            for tag, e in ev["wp"]["mass_points"].items():
                for lb, be in e["bins"].items():
                    if "threshold_loose_by_ctau" not in be:
                        raise SystemExit(f"{ev['dir']}: {tag} {lb} has no threshold_loose_by_ctau; "
                                         "rerun evaluate_workingpoint.py to add the loose WP")
                    thrs.append(_point_threshold(be)[0])
        keep_min = max(0.0, min(thrs) - 1e-4)
    print(f"[{WP_MODE}] keep-min-score {keep_min:.4g}")
    if args.rescore or not all((_wp_dir(ev) / "scores.npz").is_file() for ev in evals):
        run_scoring(evals, args, keep_min)
    for ev in evals:
        analyse(ev, args)


if __name__ == "__main__":
    main()
