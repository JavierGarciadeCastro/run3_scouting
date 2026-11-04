#!/usr/bin/env python3
import csv
import glob
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

HERE = Path(__file__).resolve().parent
FIGSIZE = (8.5, 6.5)
LXY_ORDER = ["0p0to0p2", "0p2to3p1", "3p1to11p0", "11p0to70p0"]
XSEC_H_PB = 59.8
LUMI_FB = 109.95

SRC_STYLE = {
    "table":     dict(marker="o", color="#1f77b4", mfc="#1f77b4", label="HEPData table"),
    "interp":    dict(marker="o", color="#1f77b4", mfc="white",   label="interpolated"),
    "long-ctau": dict(marker="^", color="#ff7f0e", mfc="#ff7f0e", label=r"scaled $\propto c\tau$ (beyond table)"),
    "reference": dict(marker="X", color="#7f7f7f", mfc="#7f7f7f", label="reference only (no previous limit)"),
}


def _flabel(f):
    return f"{f:g}".replace(".", "p")


def _mass_tag(r):
    if "mzd" in r:
        return f"mzd{_flabel(r['mzd'])}"
    return f"mpi{_flabel(r['mpi'])}_mA{_flabel(r['mA'])}"


def _fpr_label(entry, bin_order):
    if entry is None:
        return None
    out = []
    for lb in bin_order:
        b = entry["bins"].get(lb)
        out.append("–" if b is None else f"{-np.log10(b['fpr_target']):.0f}")
    return "/".join(out)


def make_plots(eval_dir, sub="signal_efficiency"):
    eval_dir = Path(eval_dir)
    out_dir = eval_dir / sub
    with open(out_dir / "signal_efficiency.json") as fh:
        js = json.load(fh)
    rows = js["points"]
    if not rows:
        return
    loose = js.get("wp") == "loose"
    loose_eff = js.get("loose_sig_eff", 0.9)
    wp_path = eval_dir / "majority_vote_wp.json"
    wp = json.load(open(wp_path)) if wp_path.is_file() else {"mass_points": {}}
    lxy_bins = wp.get("lxy_bins") or {}
    bin_order = [lb for lb in LXY_ORDER if lb in lxy_bins] or LXY_ORDER
    bin_ranges = " / ".join(f"{lxy_bins[lb][0]:g}–{lxy_bins[lb][1]:g}" for lb in bin_order) if lxy_bins else \
        "0–0.2 / 0.2–3.1 / 3.1–11 / 11–70"
    disc = js.get("discriminant", "BDT")
    wp_txt = f"{disc}, loose WP ({100 * loose_eff:.0f}% sig. eff.)" if loose else f"{disc}, majority-vote WP"
    for r in rows:
        r["fpr_label"] = None if loose else _fpr_label(wp["mass_points"].get(_mass_tag(r)), bin_order)

    if "mzd" in rows[0]:
        xkey, xlabel = "mzd", r"$m_{Z_d}$ [GeV]"
        groups = [(None, "HAHM")]
        def _grp(r):
            return None
    else:
        xkey, xlabel = "mpi", r"$m_{\pi_3}$ [GeV]"
        dens = sorted({int(round(r["mpi"] / r["mA"])) for r in rows})
        groups = [(d, rf"$m_{{A'}}/m_{{\pi_3}} = 1/{d}$") for d in dens]
        def _grp(r):
            return int(round(r["mpi"] / r["mA"]))
    colors = ["#1f77b4", "#d62728"]
    ctaus = sorted({r["ctau"] for r in rows})

    def _draw(ax, sel, color, label, ls="-", marker="o", annotate=None):
        sel = sorted(sel, key=lambda r: r[xkey])
        if not sel:
            return False
        ax.errorbar([r[xkey] for r in sel], [r["eff"] for r in sel], yerr=[r["err"] for r in sel],
                    color=color, ls=ls, marker=marker, markersize=5, capsize=3, linewidth=1.4, label=label)
        if annotate is not None:
            for r in sel:
                if r["fpr_label"] is None:
                    continue
                ax.annotate(r["fpr_label"], (r[xkey], r["eff"]), xytext=(0, annotate), textcoords="offset points",
                            ha="center", va="bottom" if annotate > 0 else "top", fontsize=7, color=color,
                            rotation=90, annotation_clip=False)
        return True

    def _style(ax, title):
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Signal efficiency at WP")
        ax.set_ylim(0.0, 1.05)
        ax.grid(alpha=0.3)
        ax.set_title(title, loc="left", fontsize=14)
        ax.text(1.0, 1.01, wp_txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=12)
        if loose:
            ax.axhline(loose_eff, color="gray", ls="--", lw=1, zorder=0)
        ax.tick_params(direction="in", top=True, right=True, which="both")

    def _per_ctau(sel_rows, dest, title_extra, footnote):
        dest.mkdir(parents=True, exist_ok=True)
        for c in ctaus:
            fig, ax = plt.subplots(figsize=FIGSIZE)
            drew = False
            for i, (g, lab) in enumerate(groups):
                drew |= _draw(ax, [r for r in sel_rows if r["ctau"] == c and _grp(r) == g], colors[i % 2], lab,
                              annotate=(8 if i % 2 == 0 else -8))
            if not drew:
                plt.close(fig)
                continue
            _style(ax, rf"$c\tau = {c:g}$ mm" + title_extra)
            ax.legend(loc="upper right", fontsize=12, framealpha=0.9)
            ax.text(0.01, -0.13, footnote, transform=ax.transAxes, ha="left", va="top", fontsize=9, color="dimgray")
            fig.tight_layout()
            fig.savefig(dest / f"signal_efficiency_ctau{_flabel(c)}.png", dpi=150, bbox_inches="tight")
            plt.close(fig)

    _per_ctau(rows, out_dir, "",
              rf"Score cut per $l_{{xy}}$ bin and $c\tau$ keeping {100 * loose_eff:.0f}% of the wPU-weighted test "
              rf"signal in the mass window; dashed line = {100 * loose_eff:.0f}%" if loose else
              rf"Labels: $-\log_{{10}}$(FPR target) per $l_{{xy}}$ bin [{bin_ranges}] cm")

    for lb in bin_order:
        brows = []
        for r in rows:
            d = (r.get("bins") or {}).get(lb)
            if not d or not d.get("n_in_mass_window") or d.get("eff") is None:
                continue
            frac = f"({100 * d['signal_fraction']:.0f}%)"
            brows.append({**{k: r[k] for k in ("mpi", "mA", "mzd", "ctau") if k in r},
                          "eff": d["eff"], "err": d["err"],
                          "fpr_label": frac if d.get("fpr_target") is None else
                          f"{-np.log10(d['fpr_target']):.0f} {frac}"})
        if not brows:
            continue
        rng = lxy_bins.get(lb)
        rtxt = f"{rng[0]:g} < $l_{{xy}}$ < {rng[1]:g} cm" if rng else lb
        _per_ctau(brows, out_dir / "by_lxy" / lb, f",  {rtxt}",
                  "Labels: share of the mass-window signal in this bin" if loose else
                  rf"Labels: $-\log_{{10}}$(FPR target) in this bin (share of the mass-window signal in this bin)")

    fig, ax = plt.subplots(figsize=(10, 7))
    cmap = plt.get_cmap("viridis", max(len(ctaus), 2))
    styles = [("-", "o"), ("--", "s")]
    for j, c in enumerate(ctaus):
        for i, (g, lab) in enumerate(groups):
            ls, mk = styles[i % 2]
            _draw(ax, [r for r in rows if r["ctau"] == c and _grp(r) == g], cmap(j), None, ls=ls, marker=mk)
    _style(ax, "All lifetimes")
    handles = [Line2D([], [], color=cmap(j), lw=2, label=rf"$c\tau = {c:g}$ mm") for j, c in enumerate(ctaus)]
    handles += [Line2D([], [], color="gray", ls=styles[i % 2][0], marker=styles[i % 2][1], label=lab)
                for i, (g, lab) in enumerate(groups)]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=10, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_dir / "signal_efficiency_all_ctau.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    _sr_plots(rows, groups, _grp, xkey, xlabel, bin_order, lxy_bins, wp_txt, js, out_dir)

    with open(out_dir / "signal_efficiency.json", "w") as fh:
        json.dump(js, fh, indent=1)
    print(f"Signal-efficiency plots written to {out_dir}")


def _sr_table(rows, xkey, bin_order):
    pts = {}
    for r in rows:
        key = tuple(r[k] for k in ("mpi", "mA", "mzd") if k in r)
        for lb in bin_order:
            d = (r.get("bins") or {}).get(lb)
            if d and "bkg_yield" in d:
                pts.setdefault((key, lb), {"x": r[xkey], "r": r, "d": []})["d"].append(d)
    out = {}
    for (key, lb), v in pts.items():
        ds = v["d"]
        out[(key, lb)] = {"x": v["x"], "r": v["r"],
                          "fpr_target": (float(np.median([d["fpr_target"] for d in ds]))
                                         if all(d["fpr_target"] is not None for d in ds) else None),
                          "fpr_achieved": float(np.median([d["fpr_achieved"] for d in ds])),
                          "bkg_yield": float(np.median([d["bkg_yield"] for d in ds])),
                          "bkg_raw": int(np.median([d["bkg_raw"] for d in ds])),
                          "bkg_yield_before_wp": float(np.median([d["bkg_yield_before_wp"] for d in ds]))}
    return out


def _sr_plots(rows, groups, _grp, xkey, xlabel, bin_order, lxy_bins, wp_txt, js, out_dir):
    tab = _sr_table(rows, xkey, bin_order)
    if not tab:
        print("No background yields in signal_efficiency.json (made before they were stored); "
              "rerun evaluate_workingpoint.py for the FPR / yield plots.")
        return
    lumi = js.get("lumi_fb")
    floor = js.get("b_floor", 1e-3)
    colors = ["#1f77b4", "#d62728"]
    markers = ["o", "s"]
    lumi_txt = f"{lumi:g} fb$^{{-1}}$" if lumi else ""

    def _panels(ylabel, getter, fname, annotate=None, floor_y=None, title=""):
        fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
        for lb, ax in zip(bin_order, axes.flat):
            for i, (g, lab) in enumerate(groups):
                sel = sorted((v for (k, b), v in tab.items() if b == lb and _grp(v["r"]) == g), key=lambda v: v["x"])
                if not sel:
                    continue
                xs = np.array([v["x"] for v in sel])
                ys = np.array([getter(v) for v in sel], dtype=float)
                zero = ys <= 0
                ax.plot(xs[~zero], ys[~zero], color=colors[i % 2], marker=markers[i % 2], ls="-", lw=1.4,
                        ms=6, label=lab)
                if floor_y is not None and zero.any():
                    ax.plot(xs[zero], np.full(zero.sum(), floor_y), color=colors[i % 2], marker="v", ls="none",
                            ms=8, mfc="none", label=f"{lab}: 0 MC events (floor)")
                if annotate is not None:
                    for v, y in zip(sel, np.where(zero, floor_y if floor_y else 0, ys)):
                        ax.annotate(annotate(v), (v["x"], y), xytext=(0, 7 if i % 2 == 0 else -9),
                                    textcoords="offset points", ha="center", va="bottom" if i % 2 == 0 else "top",
                                    fontsize=7, color=colors[i % 2])
            rng = lxy_bins.get(lb)
            ax.set_title(f"{rng[0]:g} < $l_{{xy}}$ < {rng[1]:g} cm" if rng else lb, loc="left", fontsize=13)
            ax.set_yscale("log")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.3, which="both")
            ax.tick_params(direction="in", top=True, right=True, which="both")
            ax.legend(fontsize=9, loc="best", framealpha=0.9)
        for ax in axes[-1]:
            ax.set_xlabel(xlabel)
        fig.suptitle(f"{wp_txt}: {title}", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
        plt.close(fig)

    if js.get("wp") == "loose":
        _panels("Background efficiency (FPR) at WP", lambda v: v["fpr_achieved"], "sr_fpr.png",
                title="MinBias FPR in mass window and lxy bin (median over ctau)")
    else:
        _panels("Background efficiency (FPR) at WP", lambda v: v["fpr_achieved"], "sr_fpr.png",
                annotate=lambda v: f"{v['fpr_target']:.0e}", title="FPR per signal region (labels: FPR target)")
    _panels(f"Expected background in SR ({lumi_txt})", lambda v: v["bkg_yield"], "sr_bkg_yield.png",
            annotate=lambda v: f"{v['bkg_raw']} raw", floor_y=floor,
            title=f"MinBias background in mass window and lxy bin after the WP, {lumi_txt} (labels: raw MC events)")
    _panels(f"Expected background in window ({lumi_txt})", lambda v: v["bkg_yield_before_wp"],
            "sr_bkg_yield_before_wp.png", title=f"MinBias background in mass window and lxy bin before the WP, {lumi_txt}")


def read_limits(path):
    out, ct = {}, None
    if not os.path.isfile(path):
        return out
    for line in open(path):
        m = re.match(r"\s*ctau=([\d.e+-]+) mm(: expected 95% CL limit r < (\S+)|:)", line)
        if m:
            ct = float(m[1])
            if m[2] != ":":
                out[ct] = float(m[3]) if m[3] != "nan" else math.nan
        m2 = re.match(r"\s*-> expected 95% CL limit r < (\S+)", line)
        if m2 and ct is not None:
            out[ct] = float(m2[1]) if m2[1] != "nan" else math.nan
    return out


def disc_label(eval_dir):
    mans = sorted(glob.glob(str(Path(eval_dir).parent / "models" / "binned_lxy" / "*_manifest.json")))
    if not mans:
        return ""
    m = json.load(open(mans[0]))
    return {("bdt", False): "BDT", ("bdt", True): "pBDT", ("pnn", False): "NN", ("pnn", True): "PNN"}[
        (m.get("model_type", "bdt"), bool(m.get("conditional", False)))]


def read_sources(eval_dir):
    rows = list(csv.DictReader(open(Path(eval_dir) / "sig_xsec_used.csv")))
    hahm = "mzd" in rows[0]
    by_tag = {}
    for r in rows:
        tag = f"mzd{_flabel(float(r['mzd']))}" if hahm else f"mpi{_flabel(float(r['mpi']))}_mA{_flabel(float(r['mA']))}"
        br = float(r["B_prev_obs"]) if r["B_prev_obs"] else None
        by_tag.setdefault(tag, {})[float(r["ctau_mm"])] = (br, r["source"])
    return hahm, by_tag


def mass_title(tag, hahm):
    if hahm:
        return "HAHM $H\\to Z_{D}Z_{D}$\n" + rf"$m_{{Z_D}} = {tag[3:].replace('p', '.')}$ GeV"
    mpi, mA = tag.split("_")
    return ("DQCD ScenarioA\n" + rf"$m_{{\pi}} = {mpi[3:].replace('p', '.')}$ GeV, "
            + rf"$m_{{A'}} = {mA[2:].replace('p', '.')}$ GeV")


def plot_one(eval_dir, tag, hahm, srcs, label):
    lim_dir = Path(eval_dir) / "best_fpr_limits"
    maj = read_limits(lim_dir / f"majority_vote_fpr_{tag}.txt")
    ctaus = sorted(set(srcs) | set(maj))
    if not ctaus:
        return None
    sigma = {c: (srcs[c][0] * XSEC_H_PB if c in srcs and srcs[c][0] else 1.0) for c in ctaus}
    abs_1pb = {c: maj[c] * sigma[c] for c in maj}

    fig, (ax, axb) = plt.subplots(2, 1, figsize=(8.5, 9.75), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.4], "hspace": 0.06})
    rs = []
    curves = ((maj, "signal = previous limit: $r$", "#1f77b4", "s"),
              (abs_1pb, r"signal = 1 pb: $\sigma\mathcal{B}$ limit [pb]", "#2ca02c", "D"))
    for d, lab, col, mk in curves:
        y = [d.get(c, math.nan) for c in ctaus]
        rs += [v for v in y if v == v and v > 0]
        ax.plot(ctaus, y, color=col, marker=mk, markersize=7, linewidth=1.9, label=lab, zorder=4)
    for c in ctaus:
        vals = [(d.get(c, math.nan), col) for d, _, col, _ in curves]
        vals = [(v, col) for v, col in vals if v == v and v > 0]
        if len(vals) == 2 and f"{vals[0][0]:.2g}" == f"{vals[1][0]:.2g}":
            vals = vals[:1]
        vals.sort(key=lambda t: -t[0])
        for i, (v, col) in enumerate(vals):
            ax.annotate(f"{v:.2g}", (c, v), textcoords="offset points", xytext=(0, 8 if i == 0 else -9),
                        ha="center", va=("bottom" if i == 0 else "top"), fontsize=8, color=col, zorder=5)
    ax.axhline(1.0, color="black", linestyle=":", linewidth=1.2, zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    if rs:
        ax.set_ylim(min(rs + [1.0]) * 0.4, max(rs + [1.0]) * 12)
    ax.set_ylabel("Expected 95% CL limit (majority-vote WP)")
    ax.set_title(label, loc="left", fontsize=15)
    ax.text(1.0, 1.01, rf"{LUMI_FB:g} fb$^{{-1}}$ (13.6 TeV, 2024)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=15)
    ax.text(0.04, 0.96, mass_title(tag, hahm), transform=ax.transAxes, va="top", ha="left", fontsize=14)
    ax.legend(loc="upper right", fontsize=12, framealpha=0.9)
    ax.tick_params(direction="in", top=True, right=True, which="both")

    have = [(c, srcs[c][0], srcs[c][1]) for c in ctaus if c in srcs and srcs[c][0]]
    if have:
        axb.plot([c for c, _, _ in have], [b for _, b, _ in have], color="0.6", linewidth=1.0, zorder=1)
    used = set()
    for c, b, src in have:
        st = SRC_STYLE.get(src, SRC_STYLE["reference"])
        axb.plot([c], [b], linestyle="none", markersize=8, markeredgewidth=1.5, zorder=3,
                 marker=st["marker"], color=st["color"], mfc=st["mfc"],
                 label=(st["label"] if src not in used else None))
        used.add(src)
        axb.annotate(f"{b:.1e}", (c, b), textcoords="offset points", xytext=(0, 7),
                     ha="center", fontsize=7, color=st["color"])
    nopb = [c for c in ctaus if c not in srcs or not srcs[c][0]]
    if nopb:
        tr = blended_transform_factory(axb.transData, axb.transAxes)
        axb.plot(nopb, [0.08] * len(nopb), transform=tr, linestyle="none", marker="v",
                 color="#d62728", markersize=8, label="not normalised (1 pb)")
    axb.set_xscale("log")
    axb.set_yscale("log")
    if have:
        lo, hi = min(b for _, b, _ in have), max(b for _, b, _ in have)
        axb.set_ylim(lo / 3.0, hi * 4.0)
    axb.set_xlabel(r"$c\tau$ [mm]")
    axb.set_ylabel((r"$\mathcal{B}_{prev}(H\to Z_{D}Z_{D})$" if hahm
                    else r"$\mathcal{B}_{prev}(H\to\psi\bar{\psi})$") + "\nused (obs.)", fontsize=12)
    sec = axb.secondary_yaxis("right", functions=(lambda b: b * XSEC_H_PB, lambda x: x / XSEC_H_PB))
    sec.set_ylabel(r"$\sigma_{sig}$ [pb]", fontsize=12)
    axb.legend(loc="best", fontsize=9, framealpha=0.9)
    axb.tick_params(direction="in", top=True, which="both")

    out_dir = Path(eval_dir) / tag.split("_")[0]
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"limit_with_prevBR_vs_ctau_{tag}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    dirs = sys.argv[1:] or sorted(glob.glob(str(HERE / "*_output" / "*" / "*" / "eval_testonly*")))
    for d in dirs:
        for sub in ("signal_efficiency", "signal_efficiency_loose"):
            if (Path(d) / sub / "signal_efficiency.json").is_file():
                make_plots(d, sub)
        if not (Path(d) / "sig_xsec_used.csv").is_file():
            print(f"skip limits of {d}: no sig_xsec_used.csv (evaluated without --sig-xsec-hepdata)")
            continue
        hahm, by_tag = read_sources(d)
        label = disc_label(d)
        n = sum(1 for tag, srcs in by_tag.items() if plot_one(d, tag, hahm, srcs, label))
        print(f"{d}: {n} limit plots")


if __name__ == "__main__":
    main()
