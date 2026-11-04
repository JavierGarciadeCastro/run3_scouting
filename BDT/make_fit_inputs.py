#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import uproot
import xgboost

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent

MFIT_MIN, MFIT_MAX = 0.1, 140.0
RAWMASS_NBINS, RAWMASS_LO, RAWMASS_HI = 15000, 0.0, 150.0

XSEC_FB = 1000.0
UNBLIND_FRAC = 0.1
ROOWEIGHT = 1.0

YEAR = "2024"
ERA = "2024"
SCENARIO = "A"
JOBS = 8
THREADS = 4

NGEN_CACHE = _HERE / "sig_ngen_cache.json"
FILTER_EFF_TABLE = _REPO / "data" / "DQCD-request.csv"

SYST_VARIANTS = ("trg_up", "trg_down", "sel_up", "sel_down")

SIG_DIR = Path("/ceph/cms/store/group/Run3Scouting/tuples_DQCD_ScenarioA")
ERA_LUMI_FB = {"C": 7.24, "D": 7.96, "E": 11.32, "H": 5.44, "I": 11.47}
LUMI_FB = round(sum(ERA_LUMI_FB.values()), 2)

_CACHE_RE = re.compile(r"Scenario(\w+?)_2024_mpi-([0-9p]+)_mA-([0-9p]+)_ctau-([0-9p]+)mm")
_CSV_RE = re.compile(r"Scenario(\w+?)_Mpi-([0-9p]+)_MA-([0-9p]+)_ctau-([0-9p]+)mm")
_SIG_RE = re.compile(r"tuples_Signal_ScenarioA_2024_mpi-([0-9p]+)_mA-([0-9p]+)_ctau-([0-9p]+)mm_2024\.root$")

def add_dxy_lxy(df):
    for sv in ("SV1", "SV2"):
        denom = df[f"{sv}_lxy"] * df[f"{sv}_mass"] / df[f"{sv}_ptmm"]
        denom = np.where(denom > 1e-9, denom, np.float32(1e-9))
        for mu in ("mu1", "mu2"):
            df[f"{sv}_{mu}_dxy_lxy"] = np.abs(df[f"{sv}_{mu}_dxy"]) / denom

def _find_binned_dir(model_arg):
    p = Path(model_arg)
    if not p.exists():
        raise SystemExit(f"model path does not exist: {p}")
    if p.is_file():
        return p.parent
    for cand in (p, p / "binned_lxy", p / "models" / "binned_lxy", p / "models"):
        if cand.is_dir() and any(cand.glob("*_manifest.json")):
            return cand
    raise SystemExit(
        f"No per-lxy models found under '{p}'. Expected '<stub>.json' + "
        f"'<stub>_manifest.json' pairs in one of: {p}, {p/'binned_lxy'}, "
        f"{p/'models'/'binned_lxy'}. Point model at the results dir written by "
        f"workingpoint.py (e.g. BDT/significance_plots_priv/ScenarioA).")

def _load_bin_models(model_arg):
    d = _find_binned_dir(model_arg)
    pairs = []
    for m in sorted(d.glob("*_manifest.json")):
        stem = m.name[:-len("_manifest.json")]
        cand = [d / (stem + ext) for ext in (".json", ".pt")]
        pairs += [c for c in cand if c.exists()][:1]
    if not pairs:
        raise SystemExit(f"No model+manifest pairs found for model {model_arg}")

    bins = []
    for mp in pairs:
        man_path = mp.parent / f"{mp.stem}_manifest.json"
        if not man_path.exists():
            raise SystemExit(f"Manifest not found next to model: {man_path}")
        with open(man_path) as f:
            man = json.load(f)
        rng = man.get("lxy_range_cm")
        if rng is None or "lxy_bin" not in man:
            raise SystemExit(
                f"{man_path.name} has no lxy_bin/lxy_range_cm -- it is not a per-lxy "
                f"model. Re-run workingpoint.py to regenerate the models.")
        bins.append({"label": man["lxy_bin"], "lo": float(rng[0]), "hi": float(rng[1]),
                     "model_path": mp, "manifest": man, "manifest_path": man_path})
    bins.sort(key=lambda b: b["lo"])

    shared = ["features", "require_l1", "resonance_veto", "hit_preselection", "conditional", "cond_vars",
              "mass_window_rel"]
    common = {}
    for k in shared:
        vals = [b["manifest"].get(k) for b in bins]
        first = vals[0]
        for b, v in zip(bins, vals):
            if v != first:
                raise SystemExit(
                    f"Per-lxy manifests disagree on '{k}': bin {bins[0]['label']} has "
                    f"{first!r} but bin {b['label']} has {v!r}. Point model at a single "
                    f"training run's output.")
        common[k] = first
    for k in ("signal_mA", "signal_mzd", "signal_ctau_mm", "signal_mratio", "scenario", "model"):
        for b in bins:
            if b["manifest"].get(k):
                common[k] = b["manifest"][k]
                break
    return bins, common

def _load_eval_wp(path):
    p = Path(path)
    if p.is_dir():
        p = p / "majority_vote_wp.json"
    if not p.exists():
        raise SystemExit(f"eval WP file not found: {p}")
    with open(p) as fh:
        return json.load(fh), p

def _eval_wp_tag(point):
    return "mpi%s_mA%s" % (_f2p(point[0]), _f2p(point[1])), point[2]

def _assign_thresholds_eval(bins, eval_wp, wp_name, point):
    wp, src = _load_eval_wp(eval_wp)
    tag, ctau = _eval_wp_tag(point)
    entry = wp.get("mass_points", {}).get(tag)
    if entry is None:
        raise SystemExit(f"No mass point '{tag}' in {src}. Available: {', '.join(sorted(wp.get('mass_points', {})))}")
    ckey = f"{ctau:g}"
    for b in bins:
        e = entry["bins"].get(b["label"])
        if e is None:
            b["threshold"], b["thr_desc"] = float("inf"), "no WP for this bin -- excluded"
            continue
        if e.get("model") and e["model"] != b["model_path"].name:
            raise SystemExit(f"{src} was made with model {e['model']} for bin {b['label']}, "
                             f"but the model points at {b['model_path'].name}.")
        if wp_name == "loose":
            if "threshold_loose_by_ctau" not in e:
                raise SystemExit(f"{src} has no threshold_loose_by_ctau for {tag} {b['label']}; "
                                 "rerun evaluate_workingpoint.py to add the loose WP.")
            thr = e["threshold_loose_by_ctau"].get(ckey)
            if thr is None:
                raise SystemExit(f"No loose threshold for ctau={ckey} mm, bin {b['label']}, point {tag} in {src}.")
            b["threshold"] = float(thr)
            b["thr_desc"] = f"loose WP sig eff={wp.get('loose_sig_eff', 0.9):g} ({src.parent.name})"
            continue
        if "threshold" in e:
            thr = e["threshold"]
        else:
            thr = e.get("threshold_by_ctau", {}).get(ckey)
            if thr is None:
                raise SystemExit(f"No threshold for ctau={ckey} mm, bin {b['label']}, point {tag} in {src}.")
        b["threshold"] = float(thr)
        b["thr_desc"] = f"majority-vote WP fpr={e['fpr_target']:g} ({src.parent.name})"
    theta = entry.get("cond_theta", {}).get(ckey)
    return wp, theta

def _p2f(s):
    return float(s.replace("p", "."))

def _f2p(v):
    return ("%g" % v).replace(".", "p")

def _same(a, b, tol=1e-3):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))

def sample_tag(scenario, mpi, mA, ctau):
    return "Signal_Scenario%s_Mpi-%s_MA-%s_ctau-%smm" % (scenario, _f2p(mpi), _f2p(mA), _f2p(ctau))

def fit_filename(fit_dir, tag):
    return str(Path(fit_dir) / ("histograms_%s_%s_%s_0.root" % (tag, ERA, YEAR)))

def lookup_ngen(scenario, mpi, mA, ctau, cache_path=NGEN_CACHE):
    cache = json.load(open(cache_path))
    hits = []
    for key in sorted(cache):
        m = _CACHE_RE.search(key)
        if not m or m.group(1) != scenario:
            continue
        if _same(_p2f(m.group(2)), mpi) and _same(_p2f(m.group(3)), mA) and _same(_p2f(m.group(4)), ctau):
            hits.append(key)
    if not hits:
        raise SystemExit("no N_gen in %s for Scenario%s mpi=%g mA=%g ctau=%g"
                         % (cache_path, scenario, mpi, mA, ctau))
    if len(hits) > 1:
        print("N_gen: %d cache entries match, using %s (others: %s)" % (len(hits), hits[0], hits[1:]))
    return float(cache[hits[0]]), hits[0]

def lookup_filter_eff(scenario, mpi, mA, ctau, table_path=FILTER_EFF_TABLE):
    exact, same_mass = [], []
    with open(table_path) as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            m = _CSV_RE.search(row[0])
            if not m or m.group(1) != scenario:
                continue
            if not (_same(_p2f(m.group(2)), mpi) and _same(_p2f(m.group(3)), mA)):
                continue
            eff = float(row[-1])
            same_mass.append(eff)
            if _same(_p2f(m.group(4)), ctau):
                exact.append(eff)
    if exact:
        return exact[0], "exact match"
    if same_mass:
        return float(np.mean(same_mass)), "mean of %d ctau points at same (mpi,mA)" % len(same_mass)
    raise SystemExit("no filter efficiency in %s for Scenario%s mpi=%g mA=%g"
                     % (table_path, scenario, mpi, mA))

def lumi_weight(is_signal, ngen=None, filter_eff=None, unblind_frac=UNBLIND_FRAC, xsec_fb=XSEC_FB, lumi_fb=LUMI_FB):
    if not is_signal:
        return 1.0
    return unblind_frac * lumi_fb * xsec_fb * filter_eff / ngen

def selection_sf(lxy):
    sf, sfup, sfdown = 1.0, 1.0, 1.0
    if lxy > 0.0 and lxy < 0.2:
        sf, sfup, sfdown = 1.10, 1.10 + 0.01, 1.10 - 0.01
    elif lxy > 0.2 and lxy < 1.0:
        sf, sfup, sfdown = 1.05, 1.05 + 0.02, 1.05 - 0.02
    elif lxy > 1.0 and lxy < 2.4:
        sf, sfup, sfdown = 1.03, 1.03 + 0.03, 1.03 - 0.03
    elif lxy > 2.4 and lxy < 3.1:
        sf, sfup, sfdown = 0.97, 0.97 + 0.06, 0.97 - 0.06
    elif lxy > 3.1 and lxy < 7.0:
        sf, sfup, sfdown = 0.82, 0.82 + 0.1, 0.82 - 0.1
    elif lxy > 7.0:
        sf, sfup, sfdown = 1.0, 1.0 + 0.3, 1.0 - 0.3
    return sf, sfup, sfdown

def trigger_sf(subpt, lxy):
    sf, sfup, sfdown = 1.0, 1.0, 1.0
    if lxy > 0.0 and lxy < 0.2:
        if subpt < 5.0:
            sf, sfup, sfdown = 1.12, 1.12 + 0.01, 1.12 - 0.01
        elif subpt > 5. and subpt < 10.:
            sf, sfup, sfdown = 0.996, 0.996 + 0.004, 0.996 - 0.003
        elif subpt > 10. and subpt < 15.:
            sf, sfup, sfdown = 0.974, 0.974 + 0.007, 0.974 - 0.008
        elif subpt > 15.0:
            sf, sfup, sfdown = 0.915, 0.915 + 0.008, 0.915 - 0.008
    elif lxy > 0.2 and lxy < 1.0:
        if subpt < 5.0:
            sf, sfup, sfdown = 1.19, 1.19 + 0.04, 1.19 - 0.05
        elif subpt > 5. and subpt < 10.:
            sf, sfup, sfdown = 1.008, 1.008 + 0.007, 1.008 - 0.008
        elif subpt > 10. and subpt < 15.:
            sf, sfup, sfdown = 0.969, 0.969 + 0.006, 0.969 - 0.008
        elif subpt > 15.0:
            sf, sfup, sfdown = 0.96, 0.96 + 0.01, 0.96 - 0.01
    elif lxy > 1.0 and lxy < 2.4:
        if subpt < 5.0:
            sf, sfup, sfdown = 1.32, 1.32 + 0.1, 1.32 - 0.2
        elif subpt > 5. and subpt < 10.:
            sf, sfup, sfdown = 1.03, 1.03 + 0.02, 1.03 - 0.02
        elif subpt > 10. and subpt < 15.:
            sf, sfup, sfdown = 0.986, 0.986 + 0.008, 0.986 - 0.017
        elif subpt > 15.0:
            sf, sfup, sfdown = 0.98, 0.98 + 0.01, 0.98 - 0.04
    elif lxy > 2.4 and subpt < 5.0:
        sf, sfup, sfdown = 1.00, 1.00 + 0.3, 1.00 - 0.3
    elif lxy > 2.4 and lxy < 3.1:
        sf, sfup, sfdown = 1.04, 1.04 + 0.01, 1.04 - 0.02
    elif lxy > 3.1 and lxy < 7.0:
        sf, sfup, sfdown = 1.02, 1.02 + 0.01, 1.02 - 0.03
    elif lxy > 7.0:
        sf, sfup, sfdown = 1.00, 1.00 + 0.05, 1.00 - 0.05
    return sf, sfup, sfdown

def per_event_mc_weight(ev):
    return ev.get("wPU", 1.0)

def _event_roow(is_signal, ev, lw):
    sf_sel = sf_sel_up = sf_sel_dn = 1.0
    sf_trg = sf_trg_up = sf_trg_dn = 1.0
    mcw = 1.0
    if is_signal:
        mcw = per_event_mc_weight(ev)
        sf_sel, sf_sel_up, sf_sel_dn = selection_sf(ev["SV1_lxy"])
        sf_trg, sf_trg_up, sf_trg_dn = trigger_sf(ev["subpt"], ev["SV1_lxy"])
    base = lw * ROOWEIGHT * mcw
    return {
        "nom":      base * sf_trg * sf_sel,
        "trg_up":   base * sf_trg_up * sf_sel,
        "trg_down": base * sf_trg_dn * sf_sel,
        "sel_up":   base * sf_trg * sf_sel_up,
        "sel_down": base * sf_trg * sf_sel_dn,
    }

def _f2p2(v):
    return ("%.2f" % v).replace(".", "p")

def dataset_name(label, wp_tag="WPTight", scenario="A", mass_point=None, mass_range=None):
    name = "Dimuon_lxy%s_%s" % (label, wp_tag)
    if mass_point:
        mpi, mA = mass_point
        name += "_Scenario%s_mPi_%s_mratio_%s" % (scenario, _f2p2(mpi), _f2p2(mA / mpi))
    if mass_range:
        name += "_%s_%s" % (_f2p2(mass_range[0]), _f2p2(mass_range[1]))
    return name

def write_fit_datasets(df, bins, out_path, is_signal, lw, wp_tag="WPTight", mass_range=None, mass_point=None):
    mlo, mhi = mass_range if mass_range else (MFIT_MIN, MFIT_MAX)
    import ROOT

    mass = df["SV1_mass"].to_numpy()
    passb = df["pass_bdt"].to_numpy()
    idx = df["lxy_bin_idx"].to_numpy()
    lxy = df["SV1_lxy"].to_numpy()
    subpt = np.minimum(df["SV1_mu1_pt"].to_numpy(), df["SV1_mu2_pt"].to_numpy())
    if is_signal and "wPU" not in df.columns:
        print("[fit] WARNING: signal input has no wPU column -- no PU weight applied")
    wpu = df["wPU"].to_numpy() if (is_signal and "wPU" in df.columns) else np.ones(len(df))

    mfit = ROOT.RooRealVar("mfit", "mfit", mlo, mhi)
    roow = ROOT.RooRealVar("roow", "roow", -1e4, 1e4)
    var_vars = {}
    if is_signal:
        for v in SYST_VARIANTS:
            var_vars[v] = ROOT.RooRealVar("roow_" + v, "roow_" + v, -1e4, 1e4)

    fout = ROOT.TFile(out_path, "RECREATE")
    print("[fit] is_signal : %s   PU weight : %s" % (is_signal, is_signal and "wPU" in df.columns))
    for i, b in enumerate(bins):
        dname = dataset_name(b["label"], wp_tag, SCENARIO, mass_point, mass_range)
        ds = ROOT.RooDataSet(dname, dname, ROOT.RooArgSet(mfit, roow), ROOT.RooFit.WeightVar("roow"))
        ds_var = {}
        if is_signal:
            for v in SYST_VARIANTS:
                ds_var[v] = ROOT.RooDataSet(dname + "_" + v, dname,
                                            ROOT.RooArgSet(mfit, var_vars[v]),
                                            ROOT.RooFit.WeightVar("roow_" + v))
        h = ROOT.TH1F(dname + "_rawmass", "; m_{#mu#mu} [GeV]; Events / 0.01 GeV",
                      RAWMASS_NBINS, RAWMASS_LO, RAWMASS_HI)

        sel = (passb > 0.5) & (idx == i) & (mass > mlo) & (mass < mhi)
        rows = np.flatnonzero(sel)
        sumw = 0.0
        for r in rows:
            ev = {"SV1_mass": float(mass[r]), "SV1_lxy": float(lxy[r]), "subpt": float(subpt[r]),
                  "wPU": float(wpu[r])}
            w = _event_roow(is_signal, ev, lw)
            mfit.setVal(ev["SV1_mass"])
            roow.setVal(w["nom"])
            ds.add(ROOT.RooArgSet(mfit, roow), w["nom"])
            h.Fill(ev["SV1_mass"], w["nom"])
            sumw += w["nom"]
            if is_signal:
                for v in SYST_VARIANTS:
                    var_vars[v].setVal(w[v])
                    ds_var[v].add(ROOT.RooArgSet(mfit, var_vars[v]), w[v])
        fout.cd()
        ds.Write()
        for v in ds_var:
            ds_var[v].Write()
        h.Write()
        print("[fit]   %-30s : %6d events, sum(w) = %.4g" % (dname, len(rows), sumw))
    fout.Close()
    print("[fit] wrote -> %s" % out_path)

def _signal_files():
    out = []
    for f in sorted(SIG_DIR.glob("*.root")):
        m = _SIG_RE.search(f.name)
        if m:
            out.append((tuple(_p2f(x) for x in m.groups()), f))
    return out

def _find_signal(files, mpi, mA, ctau):
    hits = [f for (a, b, c), f in files if _same(a, mpi) and _same(b, mA) and _same(c, ctau)]
    if len(hits) != 1:
        raise SystemExit(f"{len(hits)} signal files for mpi={mpi:g} mA={mA:g} ctau={ctau:g} in {SIG_DIR}")
    return hits[0]

def _run_job(job, log):
    env = dict(os.environ, OMP_NUM_THREADS=str(THREADS), MKL_NUM_THREADS=str(THREADS),
               OPENBLAS_NUM_THREADS=str(THREADS))
    cmd = [sys.executable, str(Path(__file__).resolve()), "--job", json.dumps(job)]
    with open(log, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    return r.returncode, log

def _g(v):
    return float(f"{v:g}")

def run_all_points(eval_dir, stage, wp_name):
    if wp_name == "loose" and stage == "score":
        raise SystemExit("--stage score is the same for both WPs (the signal files keep every event); "
                         "run it once without --wp, then --stage write --wp loose")
    if stage == "write":
        print(f"signal normalized to eras {''.join(ERA_LUMI_FB)}: {LUMI_FB:.2f} fb^-1")

    ev = Path(eval_dir).resolve()
    wp = json.load(open(ev / "majority_vote_wp.json"))
    model = wp["model_dir"]
    short = "_".join(ev.parent.name.split("_")[:2])
    root = _REPO / "fitInputs" / "2024_pnn"
    base = root / f"{short}_WP_tight"
    out = root / f"{short}_WP_{wp_name}"
    data_dir = ev / ("data_at_wp" if wp_name == "tight" else "data_at_wp_loose") / "scored"
    files = _signal_files() if stage == "score" else None

    jobs = []
    for tag, e in sorted(wp["mass_points"].items()):
        mpi, mA = float(e["mpi"]), float(e["mA"])
        names = {"mass_range": [_g(x) for x in e["mass_window"]], "mass_point": [_g(mpi), _g(mA)]}
        ctaus = sorted(float(c) for c in e.get("cond_theta", {}) or e["bins"][next(iter(e["bins"]))]["threshold_by_ctau"])
        pdir = out / tag
        sdir = base / tag / "scored"
        sdir.mkdir(parents=True, exist_ok=True)
        (pdir / "logs").mkdir(parents=True, exist_ok=True)
        if stage == "write":
            with open(pdir / "ScenarioA_2024_points.txt", "w") as fh:
                for c in ctaus:
                    fh.write(f"{mpi:g} {mA:g} {c:g}\n")
            data = data_dir / f"scored_Data2024_{tag}.root"
            if not data.is_file():
                raise SystemExit(f"missing {data}: run data_at_wp.py {ev} --wp {wp_name} --write-scored first")
            jobs.append(({"fn": "write", "model": model, "scored": str(data), "fit_dir": str(pdir), "wp_name": wp_name,
                          **names}, pdir / "logs" / "write_data.log"))
        for c in ctaus:
            point = [_g(mpi), _g(mA), _g(c)]
            scored = sdir / f"sig_ctau-{_f2p(c)}mm.root"
            if stage == "score":
                sig = _find_signal(files, mpi, mA, c)
                jobs.append(({"fn": "score", "model": model, "eval_wp": str(ev), "point": point, "input_file": str(sig),
                              "save_scored": str(scored)}, pdir / "logs" / f"score_ctau-{_f2p(c)}mm.log"))
            else:
                if not scored.is_file():
                    raise SystemExit(f"missing {scored}: run --stage score first")
                jobs.append(({"fn": "write", "model": model, "scored": str(scored), "fit_dir": str(pdir),
                              "wp_name": wp_name, "point": point, "eval_wp": str(ev) if wp_name == "loose" else None,
                              **names}, pdir / "logs" / f"write_ctau-{_f2p(c)}mm.log"))

    print(f"{stage}: {len(jobs)} jobs, {JOBS} in parallel -> {out}")
    bad = []
    with ThreadPoolExecutor(JOBS) as ex:
        for rc, log in ex.map(lambda j: _run_job(*j), jobs):
            if rc:
                bad.append(log)
                print(f"FAILED ({rc}): {log}")
    print(f"done: {len(jobs) - len(bad)} ok, {len(bad)} failed")
    sys.exit(1 if bad else 0)

def write_datasets(model, scored, fit_dir, wp_name, mass_range, mass_point, point=None, eval_wp=None):
    is_signal = point is not None
    tag = sample_tag(SCENARIO, *point) if is_signal else "Data" + YEAR
    fit_path = fit_filename(fit_dir, tag)
    Path(fit_dir).mkdir(parents=True, exist_ok=True)

    lw = 1.0
    if is_signal:
        mpi, mA, ctau = point
        ngen, key = lookup_ngen(SCENARIO, mpi, mA, ctau)
        print("[fit] N_gen       : %g  (%s)" % (ngen, key))
        feff, how = lookup_filter_eff(SCENARIO, mpi, mA, ctau)
        print("[fit] filter eff  : %g  (%s)" % (feff, how))
        lw = lumi_weight(True, ngen, feff)
        print("[fit] lumi weight : %g  = %g x %g fb^-1 x %g fb x %g / %g"
              % (lw, UNBLIND_FRAC, LUMI_FB, XSEC_FB, feff, ngen))

    bins, _common = _load_bin_models(model)
    with uproot.open(scored) as f:
        df = f["tuples"].arrays(library="pd")
    need = ["SV1_mass", "SV1_lxy", "SV1_mu1_pt", "SV1_mu2_pt", "pass_bdt", "lxy_bin_idx"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise SystemExit(f"{scored} is missing {miss}")
    if eval_wp:
        if "bdt_score" not in df.columns:
            raise SystemExit(f"{scored} has no bdt_score; cannot re-apply the {wp_name} WP")
        _assign_thresholds_eval(bins, eval_wp, wp_name, tuple(point))
        thr = np.full(len(df), np.inf, np.float32)
        idx = df["lxy_bin_idx"].to_numpy()
        for i, b in enumerate(bins):
            thr[idx == i] = np.float32(b["threshold"])
            print(f"[fit]   lxy bin {b['label']}: score > {b['threshold']:.6g}   ({b['thr_desc']})")
        df["bdt_threshold"] = thr
        df["pass_bdt"] = (df["bdt_score"].to_numpy() > thr).astype("float32")
    print(f"[fit] read {len(df)} scored events from {scored} ({int(df['pass_bdt'].sum())} pass the WP)")
    if is_signal and "sig_test_frac" in df.columns and len(df):
        tf = float(df["sig_test_frac"].iloc[0])
        lw /= tf
        print("[fit] test-only signal: lumi weight / test fraction %.4f -> %g" % (tf, lw))
    wp_tag = "WPLoose" if wp_name == "loose" else "WPTight"
    write_fit_datasets(df, bins, fit_path, is_signal, lw, wp_tag, mass_range, mass_point)

def score_signal(model, eval_wp, point, input_file, save_scored):
    point = tuple(point)
    bins, common = _load_bin_models(model)
    features = common["features"]
    require_l1 = bool(common.get("require_l1", False))
    hit_presel = bool(common.get("hit_preselection", False))
    conditional = bool(common.get("conditional", False))
    cond_vars = list(common.get("cond_vars") or [])

    _wp_file, eval_theta = _assign_thresholds_eval(bins, eval_wp, "tight", point)
    if conditional and (eval_theta is None or len(eval_theta) != len(cond_vars)):
        raise SystemExit(f"{eval_wp} has no parameter values for ctau={point[-1]:g} mm "
                         f"matching the model's cond_vars {cond_vars}.")

    needed_extra = ["SV1_mass", "SV2_mass", "SV1_lxy", "SV2_lxy", "SV1_ptmm", "SV2_ptmm",
                    "SV1_mu1_pt", "SV1_mu2_pt", "passL1", "SV1_mu1_nhitsbeforesv", "SV1_mu2_nhitsbeforesv", "wPU"]
    want = set(features) | set(needed_extra)

    print(f"[apply] models   : {len(bins)} per-lxy model(s) "
          f"({', '.join(sorted({b['model_path'].suffix for b in bins}))}) from {bins[0]['model_path'].parent}")
    print(f"[apply] features : {len(features)}"
          + (f"  (conditional: {', '.join(cond_vars)} = {eval_theta})" if conditional else ""))
    print(f"[apply] require_l1: {require_l1}")
    print(f"[apply] hit preselection (N_excess): {hit_presel}")
    print(f"[apply] working points, one per lxy bin:")
    for b in bins:
        print(f"[apply]   lxy [{b['lo']:g}, {b['hi']:g}) cm  {b['label']:>12}  "
              f"AUC {b['manifest'].get('auc', float('nan')):.3f}  "
              f"score > {b['threshold']:.6g}   ({b['thr_desc']})")

    with uproot.open(input_file) as f:
        t = f["tuples"]
        load = [b for b in want if b in set(t.keys())]
        df = t.arrays(load, library="pd").astype("float32")
        df["evtn"] = t["evtn"].array(library="np").astype(np.int64)
    n_loaded = len(df)

    if require_l1:
        if "passL1" not in df.columns:
            raise SystemExit("Model was trained with require_l1 but input has no passL1 branch "
                             "(use the tuples_L1_info tuples).")
        df = df[df["passL1"] > 0.5].reset_index(drop=True)
    if hit_presel:
        nh = df["SV1_mu1_nhitsbeforesv"].to_numpy() + df["SV1_mu2_nhitsbeforesv"].to_numpy()
        lx = df["SV1_lxy"].to_numpy()
        df = df[np.where(lx < 11.0, nh <= 0, np.where(lx < 16.0, nh <= 1, nh <= 2))].reset_index(drop=True)
    n_all = len(df)
    df = df[(df["evtn"].to_numpy() % 10) >= 7].reset_index(drop=True)
    test_frac = len(df) / n_all if n_all else 1.0
    df["sig_test_frac"] = np.float32(test_frac)
    print(f"[apply] signal test split: {len(df)} of {n_all} events (fraction {test_frac:.4f})")

    add_dxy_lxy(df)
    if conditional:
        for _c, _v in zip(cond_vars, eval_theta):
            df[_c] = np.float32(_v)
    missing = [c for c in features if c not in df.columns]
    if missing:
        raise SystemExit(
            "Input tuple is missing required feature(s) after derivation: "
            f"{missing[:10]}{'...' if len(missing) > 10 else ''}")

    lxy = df["SV1_lxy"].to_numpy()
    n_ev = len(df)
    scores = np.full(n_ev, np.nan, dtype=np.float32)
    thr_col = np.full(n_ev, np.nan, dtype=np.float32)
    bin_idx = np.full(n_ev, -1.0, dtype=np.float32)

    for i, b in enumerate(bins):
        last = (i == len(bins) - 1)
        sel = (lxy >= b["lo"]) & ((lxy <= b["hi"]) if last else (lxy < b["hi"]))
        rows = np.flatnonzero(sel)
        if len(rows) == 0:
            print(f"[apply]   lxy bin {b['label']}: no events")
            continue
        if b["model_path"].suffix == ".pt":
            if str(_HERE) not in sys.path:
                sys.path.insert(0, str(_HERE))
            from train_discriminant import PNNClassifier
            net = PNNClassifier.load(str(b["model_path"]))
            scores[rows] = net.predict_proba(
                df.iloc[rows][features].to_numpy(dtype=np.float32))[:, 1].astype("float32")
        else:
            booster = xgboost.XGBClassifier()
            booster.load_model(str(b["model_path"]))
            scores[rows] = booster.predict_proba(df.iloc[rows][features])[:, 1].astype("float32")
        thr_col[rows] = np.float32(b["threshold"])
        bin_idx[rows] = float(i)
        n_b = int((scores[rows] > b["threshold"]).sum())
        print(f"[apply]   lxy bin {b['label']}: {len(rows):,} events, {n_b:,} pass")

    df["bdt_score"] = scores
    df["bdt_threshold"] = thr_col
    df["lxy_bin_idx"] = bin_idx
    df["pass_bdt"] = (df["bdt_score"] > df["bdt_threshold"]).astype("float32")
    n_unbinned = int((bin_idx < 0).sum())
    if n_unbinned:
        print(f"[apply] {n_unbinned:,} event(s) ({100.0 * n_unbinned / max(n_ev, 1):.2f}%) have "
              f"SV1_lxy outside [{bins[0]['lo']:g}, {bins[-1]['hi']:g}] cm -- unscored.")
    print(f"[apply] loaded {n_loaded} | unbinned {n_unbinned} | pass_bdt {int(df['pass_bdt'].sum())} | kept {len(df)}")

    out_path = Path(save_scored)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(str(out_path)) as fout:
        fout.mktree("tuples", {c: df[c].values.dtype for c in df.columns})
        fout["tuples"].extend({c: df[c].values for c in df.columns})
    print(f"[apply] scored tuple ({len(df)} events) -> {out_path}")

def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--job":
        job = json.loads(sys.argv[2])
        fn = {"score": score_signal, "write": write_datasets}[job.pop("fn")]
        fn(**job)
        return
    p = argparse.ArgumentParser(
        description="Score the signal and write the per-SR RooDataSets for cpp/doAll_fitDimuonMass.C, for every "
                    "mass point and ctau of an evaluate_workingpoint.py eval dir.")
    p.add_argument("eval_dir", help="eval_testonly* dir with majority_vote_wp.json (and data_at_wp/scored/ for data).")
    p.add_argument("--stage", choices=["score", "write"], required=True,
                   help="score: score the signal tuples (BDT env, needs torch). "
                        "write: write the RooDataSets for data and signal (cmsenv, needs ROOT).")
    p.add_argument("--wp", choices=["tight", "loose"], default="tight",
                   help="Working point of the RooDataSets. loose reads the data from data_at_wp_loose/ and the "
                        "signal from the --stage score files of the tight run.")
    args = p.parse_args()
    run_all_points(args.eval_dir, args.stage, args.wp)

if __name__ == "__main__":
    main()
