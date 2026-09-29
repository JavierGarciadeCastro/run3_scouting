#!/usr/bin/env python3
"""Evaluate an already-trained scouting BDT (per-lxy-bin models produced by
workingpoint.py) WITHOUT retraining, and reproduce its ROC, discriminant,
significance and limit-vs-ctau plots.

Correctness requirement: this script rebuilds the exact same dataframe
(df_global) and the exact same train/test split as workingpoint.py (evtn % 10
split, or the random_state=42 stratified one with --do-random-splitting), so it
can identify which rows were held out as TEST during training. All plots and
yields below are then computed using ONLY the test-set rows -- never the rows
the model was trained on -- to avoid the "evaluating on training data" bias.
(workingpoint.py only trains; all evaluation happens here.)

For that reconstruction to succeed, you MUST pass this script the SAME
--scenario / --conditional / --cond-var / --mratio-grid / --holdout /
--require-l1 / --mass-point / --sig-dir / --out-name / --do-random-splitting flags that were used
to train the model you are pointing it at (these determine both the model
output directory and the exact dataframe used to compute the split). A
--plot-only flag lets you restrict which mass point(s) get PLOTTED without
touching the data-loading/split (safe to use even if the model was trained
on many mass points).

Usage (e.g., evaluate a model trained with default settings):
    python3 BDT/evaluate_workingpoint.py --out-name significance_plots_priv
"""

import argparse
import glob
import json
import os
import re
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mplhep as hep
import numpy as np
import pandas as pd
import uproot
from scipy.stats import norm
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

hep.style.use("CMS")
#plt.rcParams.update({
#    "font.size": 16, "axes.labelsize": 16, "axes.titlesize": 16,
#    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 12,
#    "legend.title_fontsize": 13, "axes.linewidth": 1.0, "xtick.major.size": 5,
#    "ytick.major.size": 5, "xtick.minor.size": 3, "ytick.minor.size": 3,
#    "xtick.major.width": 0.9, "ytick.major.width": 0.9,})

# ---------------------------------------------------------------------------
# Configuration -- must match the training run being evaluated
# ---------------------------------------------------------------------------
_parser = argparse.ArgumentParser(description="Evaluate a trained scouting BDT on its test set (no retraining).")
_parser.add_argument("--scenario", choices=["A", "B1", "B2", "C"], default="A",
                     help="Must match the --scenario used for training.")
_parser.add_argument("--holdout", nargs="+", default=[], metavar="mpi:mA[:ctau]",
                     help="Must match the --holdout used for training.")
_parser.add_argument("--conditional", action=argparse.BooleanOptionalAction, default=False,
                     help="Must match the --conditional used for training.")
_parser.add_argument("--cond-var", choices=["ctau", "mratio", "mpi"], default=["mratio"],
                     nargs="+", help="Must match the --cond-var used for training.")
_parser.add_argument("--mratio-grid", type=float, nargs="+", default=[0.33, 0.10],
                     help="Must match the --mratio-grid used for training.")
_parser.add_argument("--significance", action=argparse.BooleanOptionalAction, default=True,
                     help="Run the significance/limit evaluation (in addition to ROC/discriminant).")
_parser.add_argument("--fpr-targets", type=float, nargs="+", default=None,
                     help="Target false-positive rate(s) to evaluate (need not match training's).")
_parser.add_argument("--fpr-scan", type=float, nargs=3, default=None, metavar=("LO", "HI", "N"),
                     help="Add N log-spaced FPR targets between LO and HI (inclusive).")
_parser.add_argument("--mass-window-rel", type=float, default=0.1,
                     help="SV1 dimuon mass window half-width as a fraction r of mA for the yields. "
                          "Set to 0 to disable (full mass range).")
_parser.add_argument("--sig-dir", default="tuples_DQCD_ScenarioA",
                     help="Must match the --sig-dir used for training.")
_parser.add_argument("--out-name", default="significance_plots_priv",
                     help="Base output directory of the TRAINING run (models are read from "
                          "<out-name>.../models/binned_lxy/). Eval plots are written under "
                          "<that dir>/eval_testonly/. Ignored if --model-dir is given.")
_parser.add_argument("--model-dir", default=None, metavar="PATH",
                     help="Explicit path to the training run's Scenario output folder (the one "
                          "containing models/binned_lxy/, e.g. BDT/significance_plots_priv/"
                          "ScenarioA). Overrides --out-name (no need to reconstruct the path from "
                          "--conditional/--cond-var/--require-l1/--holdout/--scenario). The "
                          "dataset-related flags (--scenario/--conditional/--cond-var/--holdout/"
                          "--require-l1/--mass-point/--sig-dir/--mratio-grid) must still match "
                          "the training run, since they control the reconstructed train/test split.")
_parser.add_argument("--require-l1", action=argparse.BooleanOptionalAction, default=False,
                     help="Must match the --require-l1 used for training.")
_parser.add_argument("--do-random-splitting", action=argparse.BooleanOptionalAction, default=False,
                     help="Must match the --do-random-splitting used for training.")
_parser.add_argument("--mass-point", default=None, metavar="mpi:mA",
                     help="Must match the --mass-point used for training (data-loading filter). "
                          "Leave unset if training did not use --mass-point.")
_parser.add_argument("--plot-only", default=None, metavar="mpi:mA",
                     help="Optional: restrict the final plots (discriminant/significance/limit) "
                          "to this single signal mass point. Does NOT affect data loading or the "
                          "train/test split -- safe to use even for a model trained on many points.")
_args = _parser.parse_args()

def _resolve_fpr_targets():
    t = set(_args.fpr_targets or [])
    if _args.fpr_scan is not None:
        lo, hi, n = _args.fpr_scan
        t |= {float(f"{v:.3g}") for v in np.logspace(np.log10(lo), np.log10(hi), int(n))}
    if not t:
        t = {1e-1, 1e-2, 1e-3, 1e-4}
    return sorted(t, reverse=True)

def _parse_holdout(specs):
    out = []
    for s in specs:
        parts = [float(p) for p in s.split(":")]
        mpi, mA = parts[0], parts[1]
        ctau = parts[2] if len(parts) == 3 else None
        out.append((mpi, mA, ctau))
    return out

SCENARIO             = _args.scenario
use_conditional      = _args.conditional
COND_VAR             = list(_args.cond_var)
_COND_COL_OF         = {"ctau": "param_ctau", "mratio": "param_mratio", "mpi": "param_mpi"}
_COND_TAG_OF         = {"ctau": "Ctau", "mratio": "Mratio", "mpi": "Mpi"}
COND_COL             = [_COND_COL_OF[v] for v in COND_VAR]
COND_TAG             = ("_cond" + "".join(_COND_TAG_OF[v] for v in COND_VAR) if use_conditional else "")
MRATIO_GRID          = np.array(sorted(set(_args.mratio_grid)), dtype=np.float64)
FPR_TARGETS          = _resolve_fpr_targets()
# Adaptive mode kicks in only when the user didn't pin the FPR grid by hand: each (mass
# point, lxy bin) then gets its own 4 FPR targets from its own background yield in the
# window (see _adaptive_fpr_targets below), instead of sharing one grid across everything.
ADAPTIVE_FPR         = (_args.fpr_targets is None and _args.fpr_scan is None)
HOLDOUT              = _parse_holdout(_args.holdout)
HOLD_TAG             = ("_holdout-" + "-".join(s.replace(":", "-").replace(".", "p") for s in _args.holdout) if _args.holdout else "")
MASS_WINDOW_REL      = _args.mass_window_rel
MASS_WINDOW_ACTIVE   = MASS_WINDOW_REL > 0.0
REQUIRE_L1           = _args.require_l1
DO_RANDOM_SPLITTING  = _args.do_random_splitting  # True: random stratified split, False: evtn-based
FINAL_MASS_POINT     = (tuple(float(x) for x in _args.mass_point.split(":"))
                         if _args.mass_point else None)
PLOT_ONLY            = (tuple(float(x) for x in _args.plot_only.split(":"))
                         if _args.plot_only else None)
SIG_SUBDIR           = _args.sig_dir
MINBIAS_SUBDIR       = "tuples_minbias"
TUPLES_SUBDIR        = SIG_SUBDIR
LUMI_FB = 109.95
PB_TO_FB = 1.0e3
_HERE      = Path(__file__).resolve().parent
TUPLES_BASE    = Path("/ceph/cms/store/group/Run3Scouting")
tuples_dir     = TUPLES_BASE / SIG_SUBDIR
minbias_dir    = TUPLES_BASE / MINBIAS_SUBDIR
MINBIAS_FILE = "tuples_MinBias_Fil-DoubleMuOS43_2024_2024.root"
MINBIAS_NGEN_BEFOREFILTER = 8.31e9
MINBIAS_NGEN_AFTERFILTER = 409318867
MINBIAS_XSEC_BEFOREFILTER = 1.051e7
MINBIAS_XSEC = MINBIAS_XSEC_BEFOREFILTER * (MINBIAS_NGEN_AFTERFILTER / MINBIAS_NGEN_BEFOREFILTER)
# Directory workingpoint.py trained into -- this is where we READ models from. Either given
# explicitly (--model-dir) or reconstructed the same way workingpoint.py builds OUT_ROOT.
if _args.model_dir is not None:
    OUT_ROOT = Path(_args.model_dir).resolve()
else:
    OUT_ROOT = (_HERE / (_args.out_name + COND_TAG + ('_L1req' if REQUIRE_L1 else '')) / f"Scenario{SCENARIO}{HOLD_TAG}")
_BINMODEL_DIR = OUT_ROOT / "models" / "binned_lxy"
# New plots go in a subdirectory so we never overwrite the training run's own outputs.
EVAL_ROOT = OUT_ROOT / "eval_testonly"
if not _BINMODEL_DIR.is_dir():
    raise SystemExit(f"No trained models found at {_BINMODEL_DIR} -- check --model-dir (or "
                      f"--out-name/--scenario/--conditional/--cond-var/--require-l1/--holdout) "
                      f"matches the training run.")
print(f"Reading trained models from {_BINMODEL_DIR}")
if ADAPTIVE_FPR:
    print("FPR targets: adaptive (no --fpr-targets/--fpr-scan given) -- each (mass point, lxy "
          "bin) picks its own 4 targets from its own background yield in the mass window; no "
          "combined (all-lxy) significance/limit plot is produced in this mode.")
else:
    print(f"FPR targets ({len(FPR_TARGETS)}): " + ", ".join(f"{f:g}" for f in FPR_TARGETS))

def _mass_window_halfwidth(mA):
    return MASS_WINDOW_REL * mA if MASS_WINDOW_ACTIVE else None

def _mratio(mpi, mA):
    return float(MRATIO_GRID[np.argmin(np.abs(MRATIO_GRID - mA / mpi))])

def _theta_of(key):
    mpi, mA, ctau = key[0], key[1], key[2]
    _val_of = {"ctau": ctau, "mratio": _mratio(mpi, mA), "mpi": mpi}
    return tuple(_val_of[v] for v in COND_VAR)

def _wp_from_roc(fpr_arr, thr_arr, target):
    fpr_arr = np.asarray(fpr_arr, dtype=float)
    thr_arr = np.asarray(thr_arr, dtype=float)
    ok = np.flatnonzero((fpr_arr > 0.0) & (thr_arr <= 1.0))
    if len(ok) == 0:
        return None
    j = ok[int(np.argmin(np.abs(fpr_arr[ok] - target)))]
    return float(thr_arr[j]), float(fpr_arr[j])

def _holdout_mask(df):
    m = np.zeros(len(df), dtype=bool)
    if not HOLDOUT:
        return m
    is_sig = (df['label'].values == 1)
    mpi_a, mA_a, ctau_a = (df['param_mpi'].values, df['param_mA'].values, df['param_ctau'].values)
    for mpi, mA, ctau in HOLDOUT:
        sel = is_sig & np.isclose(mpi_a, mpi) & np.isclose(mA_a, mA)
        if ctau is not None:
            sel = sel & np.isclose(ctau_a, ctau)
        m |= sel
    return m

def split_by_evtn(evtn, holdout=None):
    """70/30 train/test split by event number: evtn % 10 < 7 -> train, otherwise -> test.
    Rows in `holdout` (bool mask) always go to test. Returns (train, test) row indices."""
    is_train = (np.asarray(evtn) % 10) < 7
    if holdout is not None:
        is_train &= ~holdout
    return np.flatnonzero(is_train), np.flatnonzero(~is_train)

def _p2f(s):
    return float(s.replace("p", "."))

def _flabel(f):
    return f"{f:g}".replace(".", "p")

_SIG_GLOB = f"tuples_Signal_Scenario{SCENARIO}_*2024_*.root"
_SIG_RE = re.compile(rf"tuples_Signal_Scenario{SCENARIO}_(?:(?:Par|Priv)_)?2024_mpi-(\w+)_mA-(\w+)_ctau-(\w+)mm_2024(?:_\w+)?\.root")

_sig_by_point = {}
for fpath in sorted(glob.glob(str(tuples_dir / _SIG_GLOB))):
    m = _SIG_RE.search(os.path.basename(fpath))
    if not m:
        continue
    mpi_s, mA_s, ctau_s = m.groups()
    key = (_p2f(mpi_s), _p2f(mA_s), _p2f(ctau_s))
    is_merged = bool(re.search(r"_2024\.root$", os.path.basename(fpath)))
    prev = _sig_by_point.get(key)
    if prev is None or (is_merged and not prev[1]):
        _sig_by_point[key] = (fpath, is_merged)

sig_file_params = [(f, k[0], k[1], k[2]) for k, (f, _m) in sorted(_sig_by_point.items())]
if not sig_file_params:
    raise SystemExit(f"No signal files matching {_SIG_GLOB} / {_SIG_RE.pattern} in {tuples_dir}")

if FINAL_MASS_POINT is not None:
    sig_file_params = [p for p in sig_file_params
                        if np.isclose(p[1], FINAL_MASS_POINT[0]) and np.isclose(p[2], FINAL_MASS_POINT[1])]
    if not sig_file_params:
        raise SystemExit(f"--mass-point {_args.mass_point} matched no signal point on disk")

param_grid = [(p[3], p[2], p[1]) for p in sig_file_params]

with open(_HERE / "sig_ngen_cache.json") as fh:
    _ngen_by_file = json.load(fh)

# Cache is keyed by the old Par/Priv file names: fall back to a (mpi, mA, ctau) lookup
_ngen_by_point = {}
for _fname, _n in _ngen_by_file.items():
    _m = _SIG_RE.search(_fname)
    if _m:
        _ngen_by_point.setdefault(tuple(_p2f(s) for s in _m.groups()), set()).add(_n)

SIG_NGEN = {}
for fpath, mpi_val, mA_val, ctau_val in sig_file_params:
    n = _ngen_by_file.get(os.path.basename(fpath))
    if n is None:
        _cands = _ngen_by_point.get((mpi_val, mA_val, ctau_val), set())
        if len(_cands) != 1:
            print(f"{'ambiguous' if _cands else 'no'} N_gen for {os.path.basename(fpath)} "
                  f"{sorted(_cands) if _cands else ''} -- EXCLUDED from the yield table")
            continue
        n = next(iter(_cands))
    SIG_NGEN[(mpi_val, mA_val, ctau_val)] = n

# ----------------------------------------
# BDT variables (must match workingpoint.py exactly)
# ----------------------------------------
def _make_bdt_vars():
    sv_stems = [
        "chi2Ndof", "d3d_mumu_SV", "dphi_mumu_SV", "l3d", "lxy",
        "prob", "ptmm", "x", "xErr", "y", "yErr", "z", "zErr",
        "dr_mumu", "dphi_mumu", "deta_mumu", "deta_mumu_SV",
        "sindphi_lxy", "a3d_mumu",
    ]
    mu_stems = [
        "dxy", "dxysig", "dxy_lxy", "dz", "dzsig",
        "eta", "isGlobal", "isTracker", "isvtx", "maxdr",
        "mindr", "muCSCDT", "muChambs", "muHits", "nhitsbeforesv",
        "normChi2", "phi", "phiCorr", "pixHits", "pixLayers",
        "pt", "stripHits", "trkLayers", "PFIsoAll0p3", "PFRelIsoAll0p3",
    ]
    vars_ = []
    for sv in ("SV1", "SV2"):
        for s in sv_stems:
            vars_.append(f"{sv}_{s}")
        for mu in ("mu1", "mu2"):
            for s in mu_stems:
                vars_.append(f"{sv}_{mu}_{s}")
    return vars_

BDT_VARIABLES = _make_bdt_vars()
_LOAD_BRANCHES = list(dict.fromkeys(BDT_VARIABLES + ["SV1_lxy", "SV1_mass", "SV2_mass", "passL1"]))

_READ_THREADS  = max(1, int(os.environ.get("WP_READ_THREADS", "4")))
_READ_EXECUTOR = (ThreadPoolExecutor(max_workers=_READ_THREADS) if _READ_THREADS > 1 else None)

def read_flat(path, step=2_000_000):
    with uproot.open(path) as f:
        t = f['tuples']
        available = set(t.keys())
        branches  = [b for b in _LOAD_BRANCHES if b in available]
        n = t.num_entries
        out = np.empty((n, len(branches)), dtype=np.float32)
        i = 0
        for chunk in t.iterate(branches, library='np', step_size=step, decompression_executor=_READ_EXECUTOR):
            m = len(chunk[branches[0]])
            for j, b in enumerate(branches):
                out[i:i + m, j] = chunk[b]
            i += m
        df = pd.DataFrame(out, columns=branches, copy=False)
        if not DO_RANDOM_SPLITTING:
            # Only for the train/test split (not a BDT feature). Read as int64: float32 is not exact above ~16.7M
            df['evtn'] = t['evtn'].array(library='np').astype(np.int64)
    return df

def apply_l1(df, src):
    if not REQUIRE_L1:
        return df
    return df[df['passL1'] > 0.5].reset_index(drop=True)

def compute_sample_weights(df):
    y = df['label'].values
    w = np.ones(len(df), dtype=float)
    bkg_mask = (y == 0)
    w[bkg_mask] = df.loc[bkg_mask, 'xsec_weight'].values
    n_sig   = int((y == 1).sum())
    sum_bkg = float(w[bkg_mask].sum())
    if n_sig > 0 and sum_bkg > 0:
        w[y == 1] = sum_bkg / n_sig
    return w

def add_dxy_lxy(df):
    for sv in ("SV1", "SV2"):
        denom = df[f"{sv}_lxy"] * df[f"{sv}_mass"] / df[f"{sv}_ptmm"]
        denom = np.where(denom > 1e-9, denom, np.float32(1e-9))
        for mu in ("mu1", "mu2"):
            df[f"{sv}_{mu}_dxy_lxy"] = np.abs(df[f"{sv}_{mu}_dxy"]) / denom

print("Building dataframe (10-15 minutes) -- must reproduce training's df_global exactly")
sig_frames = []
for fpath, mpi_val, mA_val, ctau_val in sig_file_params:
    df = apply_l1(read_flat(fpath), Path(fpath).name)
    df['param_ctau']   = float(ctau_val)
    df['param_mA']     = float(mA_val)
    df['param_mpi']    = float(mpi_val)
    df['param_mratio'] = _mratio(float(mpi_val), float(mA_val))
    df['label']        = 1
    sig_frames.append(df)
df_sig = pd.concat(sig_frames, ignore_index=True)
del sig_frames, df

df_bkg = apply_l1(read_flat(minbias_dir / MINBIAS_FILE), MINBIAS_FILE)
df_bkg['label'] = 0
df_bkg['xsec_weight'] = MINBIAS_XSEC / MINBIAS_NGEN_AFTERFILTER

add_dxy_lxy(df_sig)
add_dxy_lxy(df_bkg)

#lxy_bins   = [0.0, 0.2, 1.0, 2.4, 3.1, 7.0, 11.0, 16.0, 70.0]
#lxy_labels = ["0p0to0p2", "0p2to1p0", "1p0to2p4", "2p4to3p1", "3p1to7p0", "7p0to11p0", "11p0to16p0", "16p0to70p0"]

lxy_bins   = [0.0, 0.2, 3.1, 11.0, 70.0]  # Match Scouting analysis
lxy_labels = ["0p0to0p2", "0p2to3p1", "3p1to11p0", "11p0to70p0"]


lxy_pretty = {lbl: rf'[{lxy_bins[i]:g}, {lxy_bins[i+1]:g}]' for i, lbl in enumerate(lxy_labels)}

df_sig['lxy_bin'] = pd.cut(df_sig['SV1_lxy'], bins=lxy_bins, labels=lxy_labels, include_lowest=True)
df_bkg['lxy_bin'] = pd.cut(df_bkg['SV1_lxy'], bins=lxy_bins, labels=lxy_labels, include_lowest=True)

available  = set(df_sig.columns) & set(df_bkg.columns)
input_vars = [v for v in BDT_VARIABLES if v in available]
cond_vars  = list(COND_COL) if use_conditional else []

os.makedirs(EVAL_ROOT, exist_ok=True)
_roc_lxy_dir = EVAL_ROOT / "ROC_by_lxy"
os.makedirs(_roc_lxy_dir, exist_ok=True)

def _mp_dir(mpi_val, mA_val, lxy_label=None):
    d = EVAL_ROOT / f'mpi{_flabel(mpi_val)}' / f'mA_{_flabel(mA_val)}'
    if lxy_label is not None:
        d = d / f'lxy_{lxy_label}'
    os.makedirs(d, exist_ok=True)
    return d

# ---------------------------------------------------------------------------
# Reproduce the exact train/test split used by workingpoint.py (same
# --do-random-splitting, same row order) -- we only ever use the TEST half below.
# ---------------------------------------------------------------------------
if use_conditional:
    # Same placeholder theta as training (random signal theta per background event). It is
    # NOT used for the evaluation: the background test set is re-scored below at every
    # signal theta, and all of it is used against each signal point.
    _rng            = np.random.default_rng(42)
    _sig_theta_grid = df_sig[COND_COL].drop_duplicates().to_numpy()
    df_bkg[COND_COL] = _sig_theta_grid[_rng.integers(0, len(_sig_theta_grid), size=len(df_bkg))]

df_global = pd.concat([df_sig, df_bkg], ignore_index=True)
del df_sig, df_bkg

y_g   = df_global['label']
w_g   = compute_sample_weights(df_global)
df_global['weight'] = w_g
_cols  = input_vars + cond_vars
_y_all = y_g.to_numpy()

_all_i = np.arange(len(df_global))
_hold  = _holdout_mask(df_global)

if _hold.any():
    _keep_i = _all_i[~_hold]
    _hold_i = _all_i[_hold]
    if DO_RANDOM_SPLITTING:
        _tr_i, _te_i = train_test_split(
            _keep_i, test_size=0.3, random_state=42, stratify=_y_all[_keep_i])
    else:
        _tr_p, _te_p = split_by_evtn(df_global['evtn'].to_numpy()[_keep_i])   # positions within _keep_i
        _tr_i, _te_i = _keep_i[_tr_p], _keep_i[_te_p]
    _te_i = np.concatenate([_te_i, _hold_i])
else:
    if DO_RANDOM_SPLITTING:
        _tr_i, _te_i = train_test_split(_all_i, test_size=0.3, random_state=42, stratify=_y_all)
    else:
        _tr_i, _te_i = split_by_evtn(df_global['evtn'].to_numpy())

# Note: workingpoint.py's --train-max-events cap only trims _tr_i (using an
# independent RNG that doesn't touch global numpy state), so it never
# affects _te_i -- we don't need to know its value to reproduce the test set.

_col_pos  = [df_global.columns.get_loc(c) for c in _cols]
y_test_a  = _y_all[_te_i]
w_test_a  = w_g[_te_i]

_lxy_all_a = df_global['lxy_bin'].to_numpy()
lxy_test   = _lxy_all_a[_te_i]
del _lxy_all_a

_IS_TEST = np.zeros(len(df_global), dtype=bool)
_IS_TEST[_te_i] = True

print(f"Test set: {len(_te_i):,} / {len(df_global):,} rows "
      f"({100.0 * len(_te_i) / len(df_global):.1f}%) -- all plots below use these only.")

# ---------------------------------------------------------------------------
# Load the already-trained per-lxy-bin models (no fitting here)
# ---------------------------------------------------------------------------
_bin_roc_test = {}
_bin_models   = {}
_bin_thr      = {}

print("Loading trained BDTs and scoring the test set")
for lxy_label in lxy_labels:
    model_path    = _BINMODEL_DIR / f"bdt_lxy_{lxy_label}.json"
    manifest_path = _BINMODEL_DIR / f"bdt_lxy_{lxy_label}_manifest.json"
    if not model_path.is_file():
        print(f"No trained model for bin {lxy_label} -- skipping (matches training's skip).")
        continue

    with open(manifest_path) as _mf:
        _manifest = json.load(_mf)
    if _manifest.get("features") != _cols:
        print(f"[WARNING] bin {lxy_label}: feature list in manifest differs from this run's "
              f"input_vars+cond_vars -- scores may be invalid.")

    bdt_bin = XGBClassifier()
    bdt_bin.load_model(str(model_path))
    _bin_models[lxy_label] = bdt_bin

    if use_conditional:
        continue   # per-(theta, bin) ROCs are built after scoring, with the background at each theta
    m_te = (lxy_test == lxy_label)
    if m_te.sum() < 10 or len(np.unique(y_test_a[m_te])) < 2:
        print(f"bin {lxy_label}: not enough test data to build a ROC (test={int(m_te.sum())})")
        continue

    y_score_bin = bdt_bin.predict_proba(df_global.iloc[_te_i[m_te], _col_pos])[:, 1]
    fpr_t, tpr_t, thr_t = roc_curve(y_test_a[m_te], y_score_bin, sample_weight=w_test_a[m_te])
    auc_t = roc_auc_score(y_test_a[m_te], y_score_bin, sample_weight=w_test_a[m_te])
    _bin_roc_test[lxy_label] = (fpr_t, tpr_t, auc_t)
    _bin_thr[lxy_label] = (fpr_t, thr_t)
    print(f"bin {lxy_label}: AUC (test) = {auc_t:.4f}")

if not _bin_models:
    raise SystemExit("[evaluate_workingpoint] no trained models found -- nothing to evaluate.")

def _style_roc(ax, *fprs, logy=False, tprs=()):
    lo = 1e-6
    mins = [f[f > 0].min() for f in fprs if np.any(f > 0)]
    if mins:
        lo = max(1e-7, min(mins) * 0.5)
    ax.set_xscale('log')
    ax.set_xlim(lo, 1.0)
    if logy:
        ylo = 1e-3
        ymins = [t[t > 0].min() for t in tprs if np.any(t > 0)]
        if ymins:
            ylo = max(1e-4, min(ymins) * 0.5)
        ax.set_yscale('log')
        ax.set_ylim(ylo, 1.02)
    else:
        ax.set_ylim(0.0, 1.02)
    ax.grid(True, which='major', alpha=0.30, linewidth=0.6)
    ax.grid(True, which='minor', alpha=0.15, linewidth=0.4)
    return lo

# ---------------------------------------------------------------------------
# Score the TEST set only (train rows are left as NaN)
# ---------------------------------------------------------------------------
print("Scoring test-set events")
_score_cols = input_vars + cond_vars
_score_pos  = [df_global.columns.get_loc(c) for c in _score_cols]
_scores     = np.full(len(df_global), np.nan, dtype=np.float32)
_CHUNK      = 2_000_000

for lxy_label, _model in _bin_models.items():
    _rows = _te_i[lxy_test == lxy_label]
    if use_conditional:
        _rows = _rows[y_test_a[lxy_test == lxy_label] == 1]   # background: scored per theta below
    if len(_rows) == 0:
        continue
    for _i in range(0, len(_rows), _CHUNK):
        _idx = _rows[_i:_i + _CHUNK]
        _scores[_idx] = _model.predict_proba(df_global.iloc[_idx, _score_pos])[:, 1]
df_global['score'] = _scores

_SCORE_A  = df_global['score'].to_numpy()
_LABEL_A  = df_global['label'].to_numpy()
_MASS_A   = df_global['SV1_mass'].to_numpy()
_CTAU_A   = df_global['param_ctau'].to_numpy()

_lxy_ser = df_global['lxy_bin']
if not isinstance(_lxy_ser.dtype, pd.CategoricalDtype):
    _lxy_ser = _lxy_ser.astype(pd.CategoricalDtype(categories=lxy_labels))
_LXY_CODE = _lxy_ser.cat.codes.to_numpy()
_CODE_OF  = {str(_c): _i for _i, _c in enumerate(_lxy_ser.cat.categories)}
del _lxy_ser

# Signal rows, TEST-ONLY, plus the per-point fraction of reconstructed events
# that landed in test (needed below to correctly rescale efficiency, since
# SIG_NGEN counts ALL generated events, not just the test-set slice of them).
_sig_i_all = np.flatnonzero(_LABEL_A == 1)
_MPI_S = df_global['param_mpi'].to_numpy()[_sig_i_all]
_MA_S  = df_global['param_mA'].to_numpy()[_sig_i_all]
_CT_S  = _CTAU_A[_sig_i_all]
_EMPTY_I = np.empty(0, dtype=np.int64)
_SIG_ROWS = {}
SIG_TEST_FRAC = {}
for _ct, _ma, _mp in param_grid:
    _k = (_mp, _ma, _ct)
    if _k in _SIG_ROWS:
        continue
    _rows_all  = _sig_i_all[(_MPI_S == _k[0]) & (_MA_S == _k[1]) & (_CT_S == _k[2])]
    _rows_test = _rows_all[_IS_TEST[_rows_all]]
    _SIG_ROWS[_k]      = _rows_test
    SIG_TEST_FRAC[_k]  = (len(_rows_test) / len(_rows_all)) if len(_rows_all) else 0.0
del _MPI_S, _MA_S, _CT_S, _sig_i_all

# Background rows, TEST-ONLY, plus the global test fraction (used to rescale
# MINBIAS_NGEN_AFTERFILTER, which likewise counts ALL generated events).
_bkg_rows_all = np.flatnonzero(_LABEL_A == 0)
_BKG_ROWS     = _bkg_rows_all[_IS_TEST[_bkg_rows_all]]
BKG_TEST_FRAC = (len(_BKG_ROWS) / len(_bkg_rows_all)) if len(_bkg_rows_all) else 0.0
_BKG_MASS = _MASS_A[_BKG_ROWS]
del _bkg_rows_all

# Conditional BDT: score the WHOLE background test set once per signal theta (cond-var
# values set to that theta), so every signal point is evaluated against all the background.
_BKG_SCORE_THETA = {}
if use_conditional:
    _thetas = sorted({_theta_of(_k) for _k in SIG_NGEN})
    _bkg_lxy = df_global['lxy_bin'].to_numpy()[_BKG_ROWS]
    print(f"Scoring the background test set at {len(_thetas)} signal theta value(s)")
    for _theta in _thetas:
        _sc = np.full(len(_BKG_ROWS), np.nan, dtype=np.float32)
        for lxy_label, _model in _bin_models.items():
            _pos = np.flatnonzero(_bkg_lxy == lxy_label)
            for _i in range(0, len(_pos), _CHUNK):
                _p = _pos[_i:_i + _CHUNK]
                _X = df_global.iloc[_BKG_ROWS[_p], _score_pos].copy()
                for _cc, _tv in zip(COND_COL, _theta):
                    _X[_cc] = _tv
                _sc[_p] = _model.predict_proba(_X)[:, 1]
        _BKG_SCORE_THETA[_theta] = _sc
    del _bkg_lxy

def _bkg_scores(theta=None):
    """Background test-set scores aligned with _BKG_ROWS (at the given theta if conditional)."""
    return _BKG_SCORE_THETA[tuple(theta)] if use_conditional else _SCORE_A[_BKG_ROWS]

# Conditional: per-(theta, lxy bin) test-set ROC -- signal points with that theta vs the
# full background test set scored at that theta. The working-point thresholds come from these.
if use_conditional:
    _bkg_w   = w_g[_BKG_ROWS]
    _bkg_lxy = df_global['lxy_bin'].to_numpy()[_BKG_ROWS]
    _sig_lxy = df_global['lxy_bin'].to_numpy()
    for _theta in _BKG_SCORE_THETA:
        _sig_rows = np.concatenate([_r for _k, _r in _SIG_ROWS.items() if _theta_of(_k) == _theta] or [_EMPTY_I])
        for lxy_label in _bin_models:
            _rs = _sig_rows[_sig_lxy[_sig_rows] == lxy_label]
            _mb = (_bkg_lxy == lxy_label)
            if len(_rs) < 5 or _mb.sum() < 5:
                print(f"bin {lxy_label}, theta={_theta}: not enough test data to build a ROC")
                continue
            _ys = np.concatenate([np.ones(len(_rs)), np.zeros(int(_mb.sum()))])
            _ss = np.concatenate([_SCORE_A[_rs], _BKG_SCORE_THETA[_theta][_mb]])
            _ws = np.concatenate([w_g[_rs], _bkg_w[_mb]])
            fpr_t, tpr_t, thr_t = roc_curve(_ys, _ss, sample_weight=_ws)
            auc_t = roc_auc_score(_ys, _ss, sample_weight=_ws)
            _bin_roc_test.setdefault(lxy_label, {})[_theta] = (fpr_t, tpr_t, auc_t)
            _bin_thr[(_theta, lxy_label)] = (fpr_t, thr_t)
            print(f"bin {lxy_label}, theta={_theta}: AUC (test) = {auc_t:.4f}")
    del _bkg_w, _bkg_lxy, _sig_lxy

def _roc_for(lxy_label, key):
    """(fpr, thr) ROC arrays giving the working point for signal point key in lxy_label."""
    return _bin_thr.get((_theta_of(key), lxy_label) if use_conditional else lxy_label)

def _thr_key(f_t, lxy_label, key):
    return (f_t, lxy_label, _theta_of(key)) if use_conditional else (f_t, lxy_label)

def _split_by_bin(rows, scores=None):
    """Scores of rows split per lxy bin (scores aligned with rows; default: _SCORE_A[rows])."""
    out = {}
    if len(rows) == 0:
        return out
    codes  = _LXY_CODE[rows]
    scores = _SCORE_A[rows] if scores is None else scores
    for _lbl, _cd in _CODE_OF.items():
        _m = (codes == _cd)
        if _m.any():
            out[_lbl] = scores[_m]
    return out

_MASS_WINDOW_CACHE = {}
def _signal_mass_window(key):
    if not MASS_WINDOW_ACTIVE:
        return None
    if key in _MASS_WINDOW_CACHE:
        return _MASS_WINDOW_CACHE[key]
    mA   = key[1]
    half = _mass_window_halfwidth(mA)
    win  = (mA - half, mA + half)
    _MASS_WINDOW_CACHE[key] = (win, mA, half)
    return _MASS_WINDOW_CACHE[key]

_SIG_SLICES = {}
_BKG_SLICES = {}
def _sig_slices(key):
    hit = _SIG_SLICES.get(key)
    if hit is None:
        rows = _SIG_ROWS.get(key, _EMPTY_I)
        n_all = len(rows)
        res = _signal_mass_window(key)
        if res is not None and n_all:
            (lo, hi), _c, _h = res
            _m   = _MASS_A[rows]
            rows = rows[(_m >= lo) & (_m <= hi)]
        hit = (n_all, len(rows), _split_by_bin(rows))
        _SIG_SLICES[key] = hit
    return hit

def _bkg_slices(key):
    res = _signal_mass_window(key)
    lo, hi = ((None, None) if res is None else res[0])
    theta = _theta_of(key) if use_conditional else None
    ck = (lo, hi, theta)
    hit = _BKG_SLICES.get(ck)
    if hit is None:
        rows, sc = _BKG_ROWS, _bkg_scores(theta)
        if res is not None:
            _m = (_BKG_MASS >= lo) & (_BKG_MASS <= hi)
            rows, sc = rows[_m], sc[_m]
        hit = (len(rows), _split_by_bin(rows, sc))
        _BKG_SLICES[ck] = hit
    return hit

def _drop_key_slices():
    _SIG_SLICES.clear()
    _BKG_SLICES.clear()

mpi_mA_groups = {}
for ctau_val, mA_val, mpi_val in param_grid:
    mpi_mA_groups.setdefault((mpi_val, mA_val), []).append(ctau_val)

########################
###### PLOTTING ########
########################
FIGSIZE     = (8.5, 6.5)
MASS_COLORS = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4", "#e377c2"]
BKG_FACE    = "#7fc7c4"
BKG_EDGE    = "#2f5f5d"
disc_bins   = np.linspace(0.0, 1.0, 51)
disc_widths = np.diff(disc_bins)
def _theta_tex(theta):
    _tex = {"ctau": r"c\tau", "mratio": r"m_{A'}/m_{\pi_3}", "mpi": r"m_{\pi_3}"}
    return ", ".join(rf"${_tex[v]}={t:g}$" for v, t in zip(COND_VAR, theta))

for lxy_label, roc_t in _bin_roc_test.items():
    # Non-conditional: one ROC per bin. Conditional: one ROC per signal theta, with the
    # full background test set scored at that theta.
    curves = roc_t if use_conditional else {None: roc_t}
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    for i, (theta, (fpr_t, tpr_t, auc_t)) in enumerate(sorted(curves.items(), key=lambda kv: (kv[0] is None, kv[0]))):
        _lab = (f'Per-bin BDT, test set (AUC = {auc_t:.3f})' if theta is None
                else f'{_theta_tex(theta)} (AUC = {auc_t:.3f})')
        ax.plot(fpr_t, tpr_t, color=MASS_COLORS[i % len(MASS_COLORS)], linewidth=2.0, label=_lab)
    _all_fpr = [c[0] for c in curves.values()]
    _all_tpr = [c[1] for c in curves.values()]
    _lo = _style_roc(ax, *_all_fpr, logy=True, tprs=tuple(_all_tpr))
    ax.plot([_lo, 1], [_lo, 1], 'k--', alpha=0.4, linewidth=1.0)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(loc='lower right', fontsize=12 if len(curves) < 4 else 9, framealpha=0.9)
    ax.text(0.02, 0.97, "Preliminary", transform=ax.transAxes, fontsize=14, fontstyle="italic", fontweight="bold", va="top", ha="left")
    ax.tick_params(direction="in", top=True, right=True, which="both")
    fig.savefig(_roc_lxy_dir / f'ROC_bylxy_{lxy_label}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

_bkg_disc_by_bin = None if use_conditional else _split_by_bin(_BKG_ROWS)
for (mpi_val, mA_val), ctau_vals in sorted(mpi_mA_groups.items()):
    _disc_sig = {c: _split_by_bin(_SIG_ROWS.get((float(mpi_val), float(mA_val), float(c)), _EMPTY_I)) for c in ctau_vals}
    if use_conditional:
        _thetas_mp = {c: _theta_of((float(mpi_val), float(mA_val), float(c))) for c in ctau_vals}
        _bkg_disc_th = {th: _split_by_bin(_BKG_ROWS, _bkg_scores(th)) for th in set(_thetas_mp.values())}
    for lxy_label in _bin_models:
        fig, ax = plt.subplots(figsize=(7.2, 5.6))
        if not use_conditional or len(_bkg_disc_th) == 1:
            _bd = (_bkg_disc_by_bin if not use_conditional else next(iter(_bkg_disc_th.values())))
            bkg_disc = _bd.get(lxy_label, np.empty(0, dtype=np.float32))
            hb, _ = np.histogram(bkg_disc[~np.isnan(bkg_disc)], bins=disc_bins)
            if hb.sum() > 0:
                hb = hb / hb.sum()
            ax.bar(disc_bins[:-1], hb, width=disc_widths, align='edge', color=BKG_FACE, edgecolor=BKG_EDGE, linewidth=0.6, label='Background (test)', zorder=1)
        else:
            # theta differs between lifetimes (e.g. --cond-var ctau): one background curve per
            # lifetime, dashed, in that lifetime's color.
            for i, ctau_val in enumerate(sorted(ctau_vals)):
                bkg_disc = _bkg_disc_th[_thetas_mp[ctau_val]].get(lxy_label, np.empty(0, dtype=np.float32))
                hb, _ = np.histogram(bkg_disc[~np.isnan(bkg_disc)], bins=disc_bins)
                if hb.sum() > 0:
                    hb = hb / hb.sum()
                ax.stairs(hb, disc_bins, color=MASS_COLORS[i % len(MASS_COLORS)], linewidth=1.2, linestyle='--',
                          label=rf'Background at $c\tau={ctau_val:g}$ mm', zorder=2)

        _drew = False
        for i, ctau_val in enumerate(sorted(ctau_vals)):
            _sig_disc = _disc_sig[ctau_val].get(lxy_label)
            if _sig_disc is None or len(_sig_disc) < 5:
                continue
            hs, _ = np.histogram(_sig_disc, bins=disc_bins)
            if hs.sum() > 0:
                hs = hs / hs.sum()
            ax.stairs(hs, disc_bins, color=MASS_COLORS[i % len(MASS_COLORS)], linewidth=1.6, label=rf'$c\tau={ctau_val:g}$ mm', zorder=3 + i)
            _drew = True
        if not _drew:
            plt.close(fig)
            continue

        ax.set_xlabel('BDT score')
        ax.set_ylabel('a.u.')
        ax.set_xlim(0.0, 1.0)
        ax.set_title(rf"$m_{{\pi_3}} = {mpi_val:g}$ GeV, $m_{{A'}} = {mA_val:g}$ GeV")
        ax.text(0.02, 0.97, "Preliminary", transform=ax.transAxes, fontsize=14, fontstyle="italic", fontweight="bold", va="top", ha="left")
        ax.legend(loc='upper center', fontsize=12, framealpha=0.9)
        ax.tick_params(direction="in", top=True, right=True, which="both")
        _fout = (_mp_dir(mpi_val, mA_val, lxy_label) / f'Disc_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}_lxy_{lxy_label}.png')
        fig.savefig(_fout, dpi=150, bbox_inches='tight')
        plt.close(fig)
del _bkg_disc_by_bin
_bkg_disc_th = None

def _counts_above_thr(scores, t):
    return int(np.count_nonzero(scores > t))

# Efficiency/yield helpers, rescaled to account for only counting the
# TEST-set slice of reconstructed events (SIG_NGEN / MINBIAS_NGEN_AFTERFILTER
# are generator-level totals over ALL events, train+test).
def _bkg_b_at(t, key, lxy_label):
    sc = _bkg_slices(key)[1].get(lxy_label)
    if sc is None or len(sc) == 0:
        return 0.0, 0
    if BKG_TEST_FRAC <= 0.0:
        return 0.0, 0
    denom = MINBIAS_NGEN_AFTERFILTER * BKG_TEST_FRAC
    n = _counts_above_thr(sc, t)
    return MINBIAS_XSEC * PB_TO_FB * LUMI_FB * (n / denom), n

def _bkg_total_in_window(key, lxy_label):
    """Total background yield (and raw MC count) in the mass window, BEFORE any BDT
    working-point cut -- same normalization as _bkg_b_at, just with no threshold applied."""
    sc = _bkg_slices(key)[1].get(lxy_label)
    if sc is None or len(sc) == 0 or BKG_TEST_FRAC <= 0.0:
        return 0.0, 0
    denom = MINBIAS_NGEN_AFTERFILTER * BKG_TEST_FRAC
    n = len(sc)
    return MINBIAS_XSEC * PB_TO_FB * LUMI_FB * (n / denom), n

def _sig_s_at(key, t, lxy_label):
    sc = _sig_slices(key)[2].get(lxy_label)
    if sc is None or len(sc) == 0:
        return None
    frac = SIG_TEST_FRAC.get(key, 0.0)
    if frac <= 0.0:
        return None
    n_raw   = _counts_above_thr(sc, t)
    sig_eff = n_raw / (SIG_NGEN[key] * frac)
    return PB_TO_FB * LUMI_FB * sig_eff, n_raw

# ---------------------------------------------------------------------------
# Adaptive FPR grid: pick 4 targets from the NORMALIZED background yield
# (xsec * lumi, same normalization as _bkg_b_at) in the signal's mass window,
# per (mA, lxy bin) -- NOT the conditional-theta slice used later for the actual
# yield. The lowest target aims at ~O(1) normalized background event; if that
# leaves 0 raw MC events, the datacard falls back to _B_FLOOR as usual.
# ---------------------------------------------------------------------------
_BKG_LXY_CODE = _LXY_CODE[_BKG_ROWS]

def _bkg_count_in_window(mA_val, lxy_label):
    code = _CODE_OF.get(lxy_label)
    if code is None:
        return 0
    mask = (_BKG_LXY_CODE == code)
    if MASS_WINDOW_ACTIVE:
        half = _mass_window_halfwidth(mA_val)
        mask &= (_BKG_MASS >= mA_val - half) & (_BKG_MASS <= mA_val + half)
    return int(np.count_nonzero(mask))

def _bkg_yield_in_window(n_raw):
    """Normalized background yield for n_raw test-set MC events (as in _bkg_b_at, no theta coverage)."""
    if BKG_TEST_FRAC <= 0.0:
        return 0.0
    return MINBIAS_XSEC * PB_TO_FB * LUMI_FB * (n_raw / (MINBIAS_NGEN_AFTERFILTER * BKG_TEST_FRAC))

def _adaptive_fpr_targets(n_bkg):
    """4 FPR targets bracketing the normalized yield n_bkg down to ~O(1) surviving
    background event: order = ceil(log10(n_bkg)) rounds UP to the enclosing power of
    ten, then the targets are 10**-order .. 10**-(order-3) (e.g. n_bkg=1.23e7 ->
    order=8 -> 1e-8, 1e-7, 1e-6, 1e-5)."""
    if n_bkg <= 0:
        return []
    order = int(np.ceil(np.log10(n_bkg)))
    targets = {min(10.0 ** -(order - i), 1.0) for i in range(4)}
    return sorted(targets, reverse=True)

##############################
#### COMPUTE SIGNIFICANCE ####
##############################
if _args.significance:
    print("Computing significance and plotting (test set only)")

    def _thr_at_for(fpr_targets):
        """Per-(f_t, lxy bin) working point from that bin's own ROC (test set)."""
        thr_at = {}
        for f_t in fpr_targets:
            for _bk, (_fpr_b, _thr_b) in _bin_thr.items():
                # _bk is lxy_label, or (theta, lxy_label) in conditional mode
                lxy_label, _tag = (_bk[1], f" theta={_bk[0]}") if use_conditional else (_bk, "")
                _wp = _wp_from_roc(_fpr_b, _thr_b, f_t)
                if _wp is None:
                    print(f"Lxy bin {lxy_label}{_tag} has no usable cut ")
                    continue
                thr_at[(f_t, lxy_label, _bk[0]) if use_conditional else (f_t, lxy_label)] = _wp
                if _wp[1] > 3.0 * f_t:
                    print(f"Lxy bin {lxy_label}{_tag} cannot reach FPR={f_t:g}; closest achievable is {_wp[1]:.2e}")
        return thr_at

    def _fpr_rows_for(fpr_targets, thr_at, lxy_only=None):
        """(f_t, achieved_fpr, thr) rows -- median across bins, or the single bin's own
        values when lxy_only pins one bin (as in adaptive mode, where each bin has its
        own fpr_targets/thr_at to begin with)."""
        bins = [lxy_only] if lxy_only is not None else list(_bin_models)
        rows = []
        for f_t in fpr_targets:
            _wps = [wp for kk, wp in thr_at.items() if kk[0] == f_t and kk[1] in bins]
            _ts = [wp[0] for wp in _wps]
            _fa = [wp[1] for wp in _wps]
            if not _ts:
                continue
            rows.append((f_t, float(np.median(_fa)), float(np.median(_ts))))
        return rows

    if not ADAPTIVE_FPR:
        fpr_targets = list(FPR_TARGETS)
        _thr_at = _thr_at_for(fpr_targets)
        fpr_rows = _fpr_rows_for(fpr_targets, _thr_at)

    _FPR_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "*", "h", "p"]
    _FPR_COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
                   "#17becf", "#e377c2", "#8c564b", "#000000", "#bcbd22",
                   "#7f7f7f", "#aec7e8"]

    def _fpr_tex(f_t):
        e = int(np.floor(np.log10(f_t)))
        m = f_t / 10.0 ** e
        if abs(m - 10.0) < 1e-9:
            e, m = e + 1, 1.0
        return rf'10^{{{e}}}' if abs(m - 1.0) < 1e-9 else rf'{m:.3g}\times 10^{{{e}}}'

    import subprocess, uuid

    _CARDS_DIR = EVAL_ROOT / "datacards"
    os.makedirs(_CARDS_DIR, exist_ok=True)
    _LIMITS_DIR = EVAL_ROOT / "best_fpr_limits"
    os.makedirs(_LIMITS_DIR, exist_ok=True)

    def _get_asymptotic_limit(bins_sb, key, f_t, lxy_label=None, lumi_unc=1.10, workdir=str(EVAL_ROOT),
                               card_tag=None):
        mpi_val, mA_val, ctau_val = key
        bin_tag = lxy_label if lxy_label is not None else "combined"
        fpr_tag = card_tag if card_tag is not None else f"FPR{_flabel(f_t)}"
        card_name = (f"datacard_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}_ctau{_flabel(ctau_val)}"
                     f"_{fpr_tag}_{bin_tag}.txt")
        card = os.path.join(_CARDS_DIR, card_name)
        tag = uuid.uuid4().hex[:8]  # only used to keep combine's own output files unique
        n = len(bins_sb)
        bin_names = [f"bin{i + 1}" for i in range(n)]

        lines = [
            f"imax {n}", "jmax 1", "kmax *",
            "---------------------------------------------",
            "bin          " + "  ".join(bin_names),
            "observation  " + "  ".join(["-1"] * n),
            "---------------------------------------------",
            "bin          " + "  ".join(f"{bn}  {bn}" for bn in bin_names),
            "process      " + "  ".join(["sig  bkg"] * n),
            "process      " + "  ".join(["0    1"] * n),
            "rate         " + "  ".join(f"{s:g}  {b:g}" for s, b in bins_sb),
            "---------------------------------------------",
            "lumi  lnN    " + "  ".join([f"{lumi_unc}  {lumi_unc}"] * n),
        ]
        with open(card, "w") as f:
            f.write("\n".join(lines) + "\n")

        print("> combine: ", bins_sb)
        out = subprocess.run(["combine", "-M", "AsymptoticLimits", os.path.relpath(card, workdir),
                               "-n", tag, "--run", "blind",
                               "--cminDefaultMinimizerStrategy", "0",
                               "--X-rtd", "MINIMIZER_freezeDisassociatedParams",
                               "--X-rtd", "MINIMIZER_multiMin_hideConstants",
                               "--X-rtd", "MINIMIZER_multiMin_maskConstraints",
                               "--X-rtd", "MINIMIZER_multiMin_maskChannels=2"],
                              cwd=workdir, check=True, capture_output=True, text=True).stdout
        print(out)
        # datacard is kept under EVAL_ROOT/datacards/ (not deleted) so it can be inspected later.
        # combine also drops a higgsCombine<tag>.*.root per call -- clean those up, we don't need them.
        for _f in glob.glob(os.path.join(workdir, f"higgsCombine{tag}.*.root")):
            os.remove(_f)

        return {m.group(1) + "%": float(m.group(2))
                for m in re.finditer(r"Expected\s+([\d.]+)%.*?r\s*<\s*([\d.]+)", out)}

    def _significance_vs_ctau_plot(cells, keys, ctaus, mpi_val, mA_val, fpr_rows, lxy_label=None):
        if not ctaus:
            return
        order = np.argsort(np.array(ctaus, dtype=float))
        x = np.array(ctaus, dtype=float)[order]
        _MANY_FPR = len(fpr_rows) > 4

        def _span(vals, digits):
            lo, hi = min(vals), max(vals)
            def _f(v, d):
                return f'{int(v):,}' if v < 10000 and float(v).is_integer() else f'{v:.{d}g}'
            return (_f(lo, digits) if lo == hi
                    else f'{_f(lo, digits - 1)}–{_f(hi, digits - 1)}')

        fig, ax = plt.subplots(figsize=FIGSIZE)
        z_max = 0.0
        fpr_handles = []
        for i, (f_t, _f_ach, _thr) in enumerate(fpr_rows):
            z = np.array([cells[(f_t, keys[j])][1] for j in range(len(keys))], dtype=float)[order]
            z = np.nan_to_num(z, nan=0.0)
            z_max = max(z_max, float(z.max()))
            _ms = 5 if _MANY_FPR else 7
            ax.plot(x, z, color=_FPR_COLORS[i % len(_FPR_COLORS)], marker=_FPR_MARKERS[i % len(_FPR_MARKERS)], markersize=_ms, linewidth=1.4 if _MANY_FPR else 1.8, zorder=3 + i)
            _b_s  = _span([cells[(f_t, k)][3] for k in keys], 3)
            _nb_s = _span([cells[(f_t, k)][5] for k in keys], 3)
            _lab  = rf'FPR$={_fpr_tex(f_t)}$'
            _lab += (f', $B$={_b_s} ({_nb_s})' if _MANY_FPR else f'\n$B$ = {_b_s}  ({_nb_s} raw)')
            fpr_handles.append(Line2D([0], [0], color=_FPR_COLORS[i % len(_FPR_COLORS)], marker=_FPR_MARKERS[i % len(_FPR_MARKERS)], markersize=_ms, label=_lab))
        ax.set_xscale('log')
        ax.set_xlabel(r'$c\tau$ [mm]')
        ax.set_ylabel(r'Asymptotic significance $Z$')
        _ncol  = 2 if _MANY_FPR else 1
        _nrows = int(np.ceil(len(fpr_handles) / _ncol))
        _head  = 1.45 + (0.10 if _MANY_FPR else 0.15) * max(0, _nrows - 1)
        ax.set_ylim(0.0, z_max * _head if z_max > 0 else 1.0)
        ax.text(1.0, 1.01, rf'{LUMI_FB:g} fb$^{{-1}}$ (13.6 TeV, 2024)', transform=ax.transAxes, ha='right', va='bottom', fontsize=15)
        txt = [rf'Scenario {SCENARIO}', rf'$m_{{\pi_3}} = {mpi_val:g}$ GeV', rf"$m_{{A'}} = {mA_val:g}$ GeV"]
        res = _signal_mass_window(keys[0])
        if res is not None:
            (lo, hi), _c, _h = res
            txt.append(rf'Mass window: $[{lo:.2f}, {hi:.2f}]$ GeV')
        if lxy_label is not None:
            txt.append(rf'$l_{{xy}} \in {lxy_label.replace("to", "-")}$ cm')
        ax.text(0.04, 0.96, '\n'.join(txt), transform=ax.transAxes, va='top', ha='left', fontsize=14)
        leg1 = ax.legend(handles=fpr_handles, loc='upper right', framealpha=0.9,
                         fontsize=(10 if _MANY_FPR else 12), ncol=_ncol,
                         labelspacing=(0.4 if _MANY_FPR else 0.7),
                         columnspacing=0.8, handlelength=2.0,
                         title=(r'$B$ = yield (raw MC evts)' if _MANY_FPR else None),
                         title_fontsize=10)
        ax.add_artist(leg1)
        ax.tick_params(direction='in', top=True, right=True, which='both')
        fig.tight_layout()
        _fout = _mp_dir(mpi_val, mA_val, lxy_label) / f'significance_vs_ctau_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}.png'
        fig.savefig(_fout, dpi=130)
        plt.close(fig)

    def _expected_limit(cell, quantile="50.0%"):
        if len(cell) < 7 or cell[6] is None:
            return float('nan')
        return float(cell[6].get(quantile, float('nan')))

    def _best_fpr_at(cells, fpr_list, k):
        """Among fpr_list, the f_t giving the lowest (best) expected limit for key k
        (nan/non-positive limits, e.g. from s<=0 or a missing combine result, are ignored).
        Returns None if no f_t in fpr_list has a usable limit for k."""
        best_f, best_r = None, np.inf
        for f_t in fpr_list:
            r = _expected_limit(cells.get((f_t, k), ()))
            if np.isfinite(r) and r > 0.0 and r < best_r:
                best_r, best_f = r, f_t
        return best_f

    def _limit_vs_ctau_plot(cells, keys, ctaus, mpi_val, mA_val, fpr_rows, lxy_label=None):
        if not ctaus:
            return
        order = np.argsort(np.array(ctaus, dtype=float))
        x = np.array(ctaus, dtype=float)[order]
        _MANY_FPR = len(fpr_rows) > 4

        _bin_tag = lxy_label if lxy_label is not None else "all-lxy"
        for f_t, _f_ach, _thr in fpr_rows:
            for k in keys:
                r = _expected_limit(cells[(f_t, k)])
                print(f"[limit] mpi={mpi_val:g} mA={mA_val:g} ctau={k[2]:g}mm "
                      f"lxy={_bin_tag} FPR={f_t:g}: expected 95% CL limit r < {r:.4g}")

        # Total background in the mass window BEFORE any working-point cut (per key: the
        # mass window depends on the signal point).
        _totbins = [lxy_label] if lxy_label is not None else list(_bin_models)
        def _bkg_window_total(key):
            b_sum, n_sum = 0.0, 0
            for _lb in _totbins:
                b, n = _bkg_total_in_window(key, _lb)
                b_sum += b
                n_sum += n
            return b_sum, n_sum
        _bkg_tot_by_key = {k: _bkg_window_total(k) for k in keys}

        def _span(vals, digits):
            lo, hi = min(vals), max(vals)
            def _f(v, d):
                return f'{int(v):,}' if v < 10000 and float(v).is_integer() else f'{v:.{d}g}'
            return (_f(lo, digits) if lo == hi
                    else f'{_f(lo, digits - 1)}–{_f(hi, digits - 1)}')

        fig, ax = plt.subplots(figsize=FIGSIZE)
        r_max = 0.0
        r_min = np.inf
        fpr_handles = []
        for i, (f_t, _f_ach, _thr) in enumerate(fpr_rows):
            r = np.array([_expected_limit(cells[(f_t, keys[j])]) for j in range(len(keys))], dtype=float)[order]
            # nan (s<=0 or b<=0, i.e. no combine result) is left as nan on purpose -- matplotlib
            # then skips those points/segments instead of drawing a misleading r=0.
            _valid = np.isfinite(r) & (r > 0.0)
            if _valid.any():
                r_max = max(r_max, float(r[_valid].max()))
                r_min = min(r_min, float(r[_valid].min()))
            _ms = 5 if _MANY_FPR else 7
            ax.plot(x, r, color=_FPR_COLORS[i % len(_FPR_COLORS)], marker=_FPR_MARKERS[i % len(_FPR_MARKERS)], markersize=_ms, linewidth=1.4 if _MANY_FPR else 1.8, zorder=3 + i)
            _b_s  = _span([cells[(f_t, k)][3] for k in keys], 3)
            _nb_s = _span([cells[(f_t, k)][5] for k in keys], 3)
            _lab  = rf'FPR$={_fpr_tex(f_t)}$'
            _lab += (f', $B$={_b_s} ({_nb_s})' if _MANY_FPR else f'\n$B$ = {_b_s}  ({_nb_s} raw)')
            fpr_handles.append(Line2D([0], [0], color=_FPR_COLORS[i % len(_FPR_COLORS)], marker=_FPR_MARKERS[i % len(_FPR_MARKERS)], markersize=_ms, label=_lab))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$c\tau$ [mm]')
        ax.set_ylabel(r'Expected asymptotic limit on $r$ (95% CL)')
        _ncol  = 2 if _MANY_FPR else 1
        _nrows = int(np.ceil(len(fpr_handles) / _ncol))
        _head  = 1.45 + (0.10 if _MANY_FPR else 0.15) * max(0, _nrows - 1)
        if r_max > 0.0 and np.isfinite(r_min):
            ax.set_ylim(r_min * 0.5, r_max * _head)
        else:
            ax.set_ylim(1e-3, 1.0)
        ax.text(1.0, 1.01, rf'{LUMI_FB:g} fb$^{{-1}}$ (13.6 TeV, 2024)', transform=ax.transAxes, ha='right', va='bottom', fontsize=15)
        txt = [rf'Scenario {SCENARIO}', rf'$m_{{\pi_3}} = {mpi_val:g}$ GeV', rf"$m_{{A'}} = {mA_val:g}$ GeV"]
        res = _signal_mass_window(keys[0])
        if res is not None:
            (lo, hi), _c, _h = res
            txt.append(rf'Mass window: $[{lo:.2f}, {hi:.2f}]$ GeV')
        if lxy_label is not None:
            txt.append(rf'$l_{{xy}} \in {lxy_label.replace("to", "-")}$ cm')
        _bt_b_s = _span([_bkg_tot_by_key[k][0] for k in keys], 3)
        _bt_n_s = _span([_bkg_tot_by_key[k][1] for k in keys], 3)
        txt.append(rf'Bkg in window (pre-WP): {_bt_b_s} ({_bt_n_s} raw)')
        ax.text(0.04, 0.96, '\n'.join(txt), transform=ax.transAxes, va='top', ha='left', fontsize=14)
        leg1 = ax.legend(handles=fpr_handles, loc='upper right', framealpha=0.9,
                         fontsize=(10 if _MANY_FPR else 12), ncol=_ncol,
                         labelspacing=(0.4 if _MANY_FPR else 0.7),
                         columnspacing=0.8, handlelength=2.0,
                         title=(r'$B$ = yield (raw MC evts)' if _MANY_FPR else None),
                         title_fontsize=10)
        ax.add_artist(leg1)
        ax.tick_params(direction='in', top=True, right=True, which='both')
        fig.tight_layout()
        _fout = _mp_dir(mpi_val, mA_val, lxy_label) / f'limit_vs_ctau_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}.png'
        fig.savefig(_fout, dpi=130)
        plt.close(fig)

    def _bkg_mass_window_plot(mpi_val, mA_val, lxy_label, fpr_rows, thr_at, window_key):
        """Background SV1 invariant-mass shape (test set) in the signal's mass window,
        normalized to unit area: no cut at all vs. after each FPR working-point cut
        (using the same per-bin thresholds as the limit/significance plots)."""
        res = _signal_mass_window(window_key)
        if res is None or not fpr_rows:
            return
        (lo, hi), _c, _h = res
        bins_for_plot = [lxy_label] if lxy_label is not None else list(_bin_models)
        bin_codes = [_CODE_OF[b] for b in bins_for_plot if b in _CODE_OF]
        if not bin_codes:
            return
        mask_bin = np.isin(_BKG_LXY_CODE, bin_codes)
        mask_win = mask_bin & (_BKG_MASS >= lo) & (_BKG_MASS <= hi)
        mass_all = _BKG_MASS[mask_win]
        if len(mass_all) < 5:
            return

        edges = np.linspace(lo, hi, 41)
        scores_bkg = _bkg_scores(_theta_of(window_key) if use_conditional else None)

        _w_evt = _bkg_yield_in_window(1)   # xsec * lumi normalization per raw test-set MC event

        fig, ax = plt.subplots(figsize=FIGSIZE)
        h0, _ = np.histogram(mass_all, bins=edges)
        h0 = h0 * _w_evt
        ax.stairs(h0, edges, color='black', linewidth=2.2, baseline=None,
                  label=f'No cut ({h0.sum():.3g} evts, {len(mass_all):,} MC)', zorder=5)
        _h_pos = [h0[h0 > 0]]

        for i, (f_t, _fa, _thr) in enumerate(fpr_rows):
            pass_mask = np.zeros(len(_BKG_ROWS), dtype=bool)
            for lb in bins_for_plot:
                thr = thr_at.get(_thr_key(f_t, lb, window_key))
                code = _CODE_OF.get(lb)
                if thr is None or code is None:
                    continue
                pass_mask |= (_BKG_LXY_CODE == code) & (scores_bkg > thr[0])
            mass_cut = _BKG_MASS[pass_mask & mask_win]
            if len(mass_cut) < 2:
                continue
            h, _ = np.histogram(mass_cut, bins=edges)
            h = h * _w_evt
            ax.stairs(h, edges, color=_FPR_COLORS[i % len(_FPR_COLORS)], linewidth=1.8, baseline=None,
                      label=rf'FPR$={_fpr_tex(f_t)}$ ({h.sum():.3g} evts, {len(mass_cut):,} MC)', zorder=3 + i)
            _h_pos.append(h[h > 0])

        ax.set_xlabel(r'$m_{SV1}$ [GeV]')
        ax.set_ylabel(rf'Events / bin ({LUMI_FB:g} fb$^{{-1}}$)')
        ax.set_xlim(lo, hi)
        ax.set_yscale('log')
        _h_pos = np.concatenate(_h_pos)
        if len(_h_pos):
            ax.set_ylim(_h_pos.min() * 0.5, _h_pos.max() * 50)   # headroom for the legend
        ax.set_title(rf"$m_{{\pi_3}} = {mpi_val:g}$ GeV, $m_{{A'}} = {mA_val:g}$ GeV")
        ax.text(0.02, 0.97, "Preliminary", transform=ax.transAxes, fontsize=14, fontstyle="italic", fontweight="bold", va="top", ha="left")
        ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
        ax.tick_params(direction='in', top=True, right=True, which='both')
        fig.tight_layout()
        _fout = _mp_dir(mpi_val, mA_val, lxy_label) / f'bkg_mass_window_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}.png'
        fig.savefig(_fout, dpi=130)
        plt.close(fig)

    def _best_fpr_heatmap(best_map, lxy_order, ctaus, mpi_val, mA_val):
        """2D map (lxy bin x ctau) of the FPR that gives the best (lowest) expected limit
        in each cell, for one signal mass point. best_map: {(lxy_label, ctau): best_f_t}."""
        ctaus_sorted = sorted(set(ctaus), key=float)
        lxy_used = [l for l in lxy_order if any((l, c) in best_map for c in ctaus_sorted)]
        if not ctaus_sorted or not lxy_used:
            return
        n_x, n_y = len(lxy_used), len(ctaus_sorted)
        grid = np.full((n_y, n_x), np.nan)
        for iy, ct in enumerate(ctaus_sorted):
            for ix, lxy_label in enumerate(lxy_used):
                f = best_map.get((lxy_label, ct))
                if f is not None:
                    grid[iy, ix] = np.log10(f)

        fig, ax = plt.subplots(figsize=(10, 10),
                                constrained_layout=True)
        im = ax.imshow(grid, aspect='auto', origin='lower', cmap='viridis_r')
        ax.set_xticks(range(n_x))
        ax.set_xticklabels([lxy_pretty[l] for l in lxy_used], rotation=40, ha='right')
        ax.set_yticks(range(n_y))
        ax.set_yticklabels([f'{c:g}' for c in ctaus_sorted])
        ax.set_xlabel(r'$l_{xy}$ bin [cm]')
        ax.set_ylabel(r'$c\tau$ [mm]')
        _finite = grid[np.isfinite(grid)]
        _mid = 0.5 * (_finite.min() + _finite.max()) if _finite.size else 0.0
        for iy in range(n_y):
            for ix in range(n_x):
                f = best_map.get((lxy_used[ix], ctaus_sorted[iy]))
                txt = f'{f:.1e}' if f is not None else '--'
                color = 'white' if (np.isfinite(grid[iy, ix]) and grid[iy, ix] < _mid) else 'black'
                ax.text(ix, iy, txt, ha='center', va='center', fontsize=11, color=color)
        cb = fig.colorbar(im, ax=ax)
        cb.set_label(r'Best FPR ($\log_{10}$ scale)')
        ax.set_title(rf"$m_{{\pi_3}} = {mpi_val:g}$ GeV, $m_{{A'}} = {mA_val:g}$ GeV")
        ax.text(0.02, 1.03, "Preliminary", transform=ax.transAxes, fontsize=13, fontstyle="italic", fontweight="bold", va="bottom", ha="left")
        _fout = _mp_dir(mpi_val, mA_val) / f'best_fpr_heatmap_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}.png'
        fig.savefig(_fout, dpi=130)
        plt.close(fig)

    _B_FLOOR = 1e-3  # B=0 is treated as this floor instead of skipping the point (s<=0 still skips)

    def _compute_significance_cells(keys, thr_at, fpr_rows, lxy_only=None):
        _bins = [lxy_only] if lxy_only is not None else list(_bin_models)
        cells = {}
        for f_t, _f_ach, _thr_display in fpr_rows:
            for k in keys:
                s_tot, b_tot, saw_s = 0.0, 0.0, False
                ns_tot, nb_tot = 0, 0
                bins_sb = []
                for _lb in _bins:
                    thr = thr_at.get(_thr_key(f_t, _lb, k))
                    if thr is None:
                        continue
                    thr = thr[0]
                    s_bin = 0.0
                    s = _sig_s_at(k, thr, _lb)
                    b, n_b = _bkg_b_at(thr, k, _lb)
                    if s is not None:
                        s_bin, n_s = s
                        s_tot += s_bin
                        ns_tot += n_s
                        saw_s = True
                    b_tot  += b
                    nb_tot += n_b
                    bins_sb.append((s_bin, b if b > 0.0 else _B_FLOOR))
                s = s_tot if saw_s else None
                b = b_tot if b_tot > 0.0 else _B_FLOOR
                if s is None or s <= 0.0:
                    cells[(f_t, k)] = (float('nan'), 0.0, s, b, ns_tot, nb_tot)
                    continue
                Z  = float(np.sqrt(2.0 * ((s + b) * np.log1p(s / b) - s)))
                p0 = float(norm.sf(Z))
                limit = _get_asymptotic_limit(bins_sb, key=k, f_t=f_t, lxy_label=lxy_only)
                cells[(f_t, k)] = (p0, Z, s, b, ns_tot, nb_tot, limit)
        return cells

    def _combined_limit_with_bin_fpr(k, fpr_for_bin, card_tag):
        """Combine all lxy bins for signal point k, each cut at its OWN f_t from
        fpr_for_bin ({lxy_label: f_t}) instead of one shared f_t -- one extra combine
        call. Returns (expected 50% limit or nan, {lxy_label: f_t actually used})."""
        bins_sb, used = [], {}
        for lxy_label in _bin_models:
            f_t = fpr_for_bin.get(lxy_label)
            if f_t is None:
                continue
            _bin_roc = _roc_for(lxy_label, k)
            if _bin_roc is None:
                continue
            fpr_b, thr_b = _bin_roc
            wp = _wp_from_roc(fpr_b, thr_b, f_t)
            if wp is None:
                continue
            thr = wp[0]
            s = _sig_s_at(k, thr, lxy_label)
            s_bin = s[0] if s is not None else 0.0
            b, _n_b = _bkg_b_at(thr, k, lxy_label)
            bins_sb.append((s_bin, b if b > 0.0 else _B_FLOOR))
            used[lxy_label] = f_t
        if not bins_sb or not any(s_bin > 0.0 for s_bin, _ in bins_sb):
            return float('nan'), used
        limit_dict = _get_asymptotic_limit(bins_sb, key=k, f_t=0.0, lxy_label="combined", card_tag=card_tag)
        return float(limit_dict.get("50.0%", float('nan'))), used

    def _best_fpr_strategy_plot(mpi_val, mA_val, keys, majority_r, perlifetime_r,
                                 fixed_cells=None, fpr_rows=None):
        """Limit vs ctau comparing the two best-FPR-per-bin re-combine strategies (and, in
        fixed-grid mode, the individual shared-FPR curves already in fixed_cells) on one plot."""
        ctaus = [k[2] for k in keys]
        if not ctaus:
            return
        order = np.argsort(np.array(ctaus, dtype=float))
        x = np.array(ctaus, dtype=float)[order]

        fig, ax = plt.subplots(figsize=FIGSIZE)
        r_max, r_min = 0.0, np.inf

        def _track(r):
            nonlocal r_max, r_min
            _valid = np.isfinite(r) & (r > 0.0)
            if _valid.any():
                r_max = max(r_max, float(r[_valid].max()))
                r_min = min(r_min, float(r[_valid].min()))

        if fixed_cells is not None and fpr_rows:
            # Offset by 2 so the fixed-grid colors never collide with the bold blue/red
            # majority-vote/per-lifetime curves plotted below.
            _grid_off = 2
            for i, (f_t, _fa, _thr) in enumerate(fpr_rows):
                r = np.array([_expected_limit(fixed_cells.get((f_t, k), ())) for k in keys], dtype=float)[order]
                _track(r)
                ax.plot(x, r, color=_FPR_COLORS[(i + _grid_off) % len(_FPR_COLORS)],
                         marker=_FPR_MARKERS[(i + _grid_off) % len(_FPR_MARKERS)],
                         markersize=5, linewidth=1.2, linestyle='--', alpha=0.85, zorder=2,
                         label=rf'FPR$={_fpr_tex(f_t)}$ (shared)')

        for r_by_ctau, label, color, marker in (
                (majority_r, 'Majority-vote FPR per bin', '#1f77b4', 's'),
                (perlifetime_r, 'Per-lifetime best FPR per bin', '#d62728', 'o')):
            r = np.array([r_by_ctau.get(k[2], float('nan')) for k in keys], dtype=float)[order]
            _track(r)
            ax.plot(x, r, color=color, marker=marker, markersize=7, linewidth=1.9, zorder=4, label=label)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$c\tau$ [mm]')
        ax.set_ylabel(r'Expected asymptotic limit on $r$ (95% CL)')
        if r_max > 0.0 and np.isfinite(r_min):
            ax.set_ylim(r_min * 0.5, r_max * 1.6)
        else:
            ax.set_ylim(1e-3, 1.0)
        ax.text(1.0, 1.01, rf'{LUMI_FB:g} fb$^{{-1}}$ (13.6 TeV, 2024)', transform=ax.transAxes,
                ha='right', va='bottom', fontsize=15)
        txt = [rf'Scenario {SCENARIO}', rf'$m_{{\pi_3}} = {mpi_val:g}$ GeV', rf"$m_{{A'}} = {mA_val:g}$ GeV"]
        ax.text(0.04, 0.96, '\n'.join(txt), transform=ax.transAxes, va='top', ha='left', fontsize=14)
        ax.legend(loc='upper right', fontsize=12, framealpha=0.9)
        ax.tick_params(direction='in', top=True, right=True, which='both')
        fig.tight_layout()
        _fout = _mp_dir(mpi_val, mA_val) / f'best_fpr_strategy_limit_vs_ctau_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}.png'
        fig.savefig(_fout, dpi=130)
        plt.close(fig)

    def _majority_vote_fpr(best_map, lxy_order, ctaus):
        """Per lxy bin, the FPR value voted most often across ctaus by best_map[(lxy,ctau)]
        (ties broken by the tightest/smallest FPR). Bins with no vote at all are omitted."""
        vote = {}
        for lxy_label in lxy_order:
            counts = Counter(best_map[(lxy_label, c)] for c in ctaus if (lxy_label, c) in best_map)
            if not counts:
                continue
            top = counts.most_common()
            best_n = top[0][1]
            winners = sorted(f for f, n in top if n == best_n)
            vote[lxy_label] = (winners[0], counts)
        return vote

    _final_groups = mpi_mA_groups
    if PLOT_ONLY is not None:
        _final_groups = {k: v for k, v in mpi_mA_groups.items()
                          if np.isclose(k[0], PLOT_ONLY[0]) and np.isclose(k[1], PLOT_ONLY[1])}
        if not _final_groups:
            print(f"--plot-only {_args.plot_only} matched no signal point; "
                  f"available: {sorted(mpi_mA_groups)}")

    # Majority-vote working points, machine-readable, for a downstream analysis to apply:
    # for mass point (mpi, mA) and lxy bin, keep events with
    #   predict_proba(bdt_lxy_<bin>.json)[:, 1] > threshold   (and SV1_mass in mass_window).
    # Merged with an existing file so --plot-only runs only update their own mass point.
    _WP_PATH = EVAL_ROOT / "majority_vote_wp.json"
    _wp_out = {}
    if _WP_PATH.is_file():
        try:
            with open(_WP_PATH) as fh:
                _wp_out = json.load(fh)
        except (OSError, json.JSONDecodeError):
            _wp_out = {}
    _wp_out.update({
        "description": "Majority-vote BDT working point per (mpi, mA) and lxy bin. Cut: "
                       "score > threshold, score = predict_proba(model)[:, 1] with 'features' as "
                       "input; events in lxy bin (lo, hi] of SV1_lxy [cm] (first bin includes lo) "
                       "and SV1_mass in mass_window [GeV]. Conditional BDT: set the cond_vars "
                       "features of EVERY event (signal and background) to cond_theta[ctau] and "
                       "cut at threshold_by_ctau[ctau].",
        "scenario": SCENARIO,
        "model_dir": str(_BINMODEL_DIR),
        "features": list(_cols),
        "conditional": bool(use_conditional),
        "cond_vars": list(COND_VAR) if use_conditional else [],
        "require_l1": bool(REQUIRE_L1),
        "mass_window_rel": MASS_WINDOW_REL if MASS_WINDOW_ACTIVE else None,
        "lxy_bins": {lb: [lxy_bins[i], lxy_bins[i + 1]] for i, lb in enumerate(lxy_labels)},
        "b_floor": _B_FLOOR,
    })
    _wp_out.setdefault("mass_points", {})

    for (mpi_val, mA_val), ctau_vals in sorted(_final_groups.items()):
        keys = [(mpi_val, mA_val, c) for c in sorted(ctau_vals) if (mpi_val, mA_val, c) in SIG_NGEN]
        if not keys:
            continue

        ctaus = [k[2] for k in keys]
        best_fpr_map = {}  # (lxy_label, ctau) -> f_t giving the lowest expected limit

        if not ADAPTIVE_FPR:
            # Fixed FPR grid, shared across all lxy bins -> a combined (all-lxy) point makes
            # sense, since every bin is being cut at the same nominal FPR.
            cells = _compute_significance_cells(keys, _thr_at, fpr_rows)
            _significance_vs_ctau_plot(cells, keys, ctaus, mpi_val, mA_val, fpr_rows)
            _limit_vs_ctau_plot(cells, keys, ctaus, mpi_val, mA_val, fpr_rows)
            _bkg_mass_window_plot(mpi_val, mA_val, None, fpr_rows, _thr_at, keys[0])

            _fpr_list = [row[0] for row in fpr_rows]
            for lxy_label in _bin_models:
                cbl = _compute_significance_cells(keys, _thr_at, fpr_rows, lxy_only=lxy_label)
                _significance_vs_ctau_plot(cbl, keys, ctaus, mpi_val, mA_val, fpr_rows, lxy_label=lxy_label)
                _limit_vs_ctau_plot(cbl, keys, ctaus, mpi_val, mA_val, fpr_rows, lxy_label=lxy_label)
                _bkg_mass_window_plot(mpi_val, mA_val, lxy_label, fpr_rows, _thr_at, keys[0])
                for k in keys:
                    best_f = _best_fpr_at(cbl, _fpr_list, k)
                    if best_f is not None:
                        best_fpr_map[(lxy_label, k[2])] = best_f
        else:
            # Adaptive: each lxy bin gets its own FPR grid from its own background yield in
            # the mass window, so bins are no longer cut at a shared nominal FPR -- there is
            # no meaningful combined (all-lxy) point in this mode, only per-bin plots.
            for lxy_label in _bin_models:
                n_bkg_raw = _bkg_count_in_window(mA_val, lxy_label)
                n_bkg = _bkg_yield_in_window(n_bkg_raw)
                bin_fpr_targets = _adaptive_fpr_targets(n_bkg)
                if not bin_fpr_targets:
                    print(f"mpi={mpi_val:g} mA={mA_val:g} lxy={lxy_label}: no background in the "
                          f"mass window (n={n_bkg_raw}) -- skipping adaptive FPR scan for this bin.")
                    continue
                thr_at_bin = _thr_at_for(bin_fpr_targets)
                fpr_rows_bin = _fpr_rows_for(bin_fpr_targets, thr_at_bin, lxy_only=lxy_label)
                if not fpr_rows_bin:
                    continue
                print(f"mpi={mpi_val:g} mA={mA_val:g} lxy={lxy_label}: B(window)={n_bkg:.3g} ({n_bkg_raw} raw) "
                      f"-> adaptive FPR targets " + ", ".join(f"{f:g}" for f in bin_fpr_targets))
                cbl = _compute_significance_cells(keys, thr_at_bin, fpr_rows_bin, lxy_only=lxy_label)
                _significance_vs_ctau_plot(cbl, keys, ctaus, mpi_val, mA_val, fpr_rows_bin, lxy_label=lxy_label)
                _limit_vs_ctau_plot(cbl, keys, ctaus, mpi_val, mA_val, fpr_rows_bin, lxy_label=lxy_label)
                _bkg_mass_window_plot(mpi_val, mA_val, lxy_label, fpr_rows_bin, thr_at_bin, keys[0])
                for k in keys:
                    best_f = _best_fpr_at(cbl, bin_fpr_targets, k)
                    if best_f is not None:
                        best_fpr_map[(lxy_label, k[2])] = best_f

        _bin_order = [l for l in lxy_labels if l in _bin_models]
        _best_fpr_heatmap(best_fpr_map, _bin_order, ctaus, mpi_val, mA_val)

        # ---- Extra combine re-runs: fix a FPR-per-bin choice, then combine again ----

        # (1) Majority vote: one FPR per bin, chosen by the mode across lifetimes of the
        # per-(bin,lifetime) best FPR, then re-combined per lifetime with that single,
        # lifetime-independent cut for each bin.
        majority_r, perlifetime_r = {}, {}

        vote = _majority_vote_fpr(best_fpr_map, _bin_order, ctaus)
        if vote:
            fpr_for_bin = {lb: f for lb, (f, _c) in vote.items()}
            _mv_path = _LIMITS_DIR / f'majority_vote_fpr_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}.txt'
            with open(_mv_path, "w") as fh:
                fh.write(f"Mass point mpi={mpi_val:g} mA={mA_val:g}\n")
                fh.write("Majority-vote FPR per lxy bin (mode of the per-(bin,lifetime) best FPR "
                         "across all lifetimes; ties broken by the tightest FPR):\n")
                for lb in _bin_order:
                    if lb not in vote:
                        fh.write(f"  {lb}: no vote (no valid best FPR for any lifetime)\n")
                        continue
                    f_win, counts = vote[lb]
                    votes_str = ", ".join(f"{f:g}:{n}" for f, n in sorted(counts.items()))
                    fh.write(f"  {lb}: FPR = {f_win:g}  (votes {votes_str})\n")
                fh.write("\nCombined limit re-run with that fixed per-bin FPR, per lifetime:\n")
                for k in keys:
                    r, _used = _combined_limit_with_bin_fpr(k, fpr_for_bin, card_tag="majvote")
                    majority_r[k[2]] = r
                    fh.write(f"  ctau={k[2]:g} mm: expected 95% CL limit r < {r:.4g}\n")
            print(f"Wrote majority-vote FPR summary to {_mv_path}")

            _wp_bins = {}
            for lb in _bin_order:
                if lb not in vote:
                    continue
                f_win, counts = vote[lb]
                _entry = {
                    "model": f"bdt_lxy_{lb}.json",
                    "fpr_target": float(f_win),
                    "votes": {f"{f:g}": int(n) for f, n in sorted(counts.items())},
                }
                if not use_conditional:
                    roc = _bin_thr.get(lb)
                    wp = _wp_from_roc(*roc, f_win) if roc is not None else None
                    if wp is None:
                        continue
                    _entry["fpr_achieved"] = wp[1]
                    _entry["threshold"] = wp[0]
                else:
                    # the threshold depends on the signal theta -> one per lifetime (see cond_theta)
                    _thr_c, _fpr_c = {}, {}
                    for k in keys:
                        roc = _roc_for(lb, k)
                        wp = _wp_from_roc(*roc, f_win) if roc is not None else None
                        if wp is not None:
                            _thr_c[f"{k[2]:g}"], _fpr_c[f"{k[2]:g}"] = wp[0], wp[1]
                    if not _thr_c:
                        continue
                    _entry["fpr_achieved_by_ctau"] = _fpr_c
                    _entry["threshold_by_ctau"] = _thr_c
                _wp_bins[lb] = _entry
            _res = _signal_mass_window(keys[0])
            _wp_entry = {
                "mpi": float(mpi_val),
                "mA": float(mA_val),
                "mass_window": [float(x) for x in _res[0]] if _res is not None else None,
                "bins": _wp_bins,
                "expected_limit_r": {f"{c:g}": (float(r) if np.isfinite(r) else None)
                                     for c, r in majority_r.items()},
            }
            if use_conditional:
                _wp_entry["cond_theta"] = {f"{k[2]:g}": [float(v) for v in _theta_of(k)] for k in keys}
            _wp_out["mass_points"][f"mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}"] = _wp_entry
            with open(_WP_PATH, "w") as fh:
                json.dump(_wp_out, fh, indent=2)
            print(f"Wrote majority-vote working points to {_WP_PATH}")

        # (2) Per-lifetime best: each bin uses its own best FPR for that specific lifetime
        # (may differ bin-to-bin and lifetime-to-lifetime).
        _pl_path = _LIMITS_DIR / f'per_lifetime_best_fpr_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}.txt'
        with open(_pl_path, "w") as fh:
            fh.write(f"Mass point mpi={mpi_val:g} mA={mA_val:g}\n")
            fh.write("Per-lifetime combined limit re-run: each bin cut at its OWN best FPR for "
                     "that specific lifetime:\n")
            for k in keys:
                fpr_for_bin = {lb: best_fpr_map[(lb, k[2])] for lb in _bin_order if (lb, k[2]) in best_fpr_map}
                r, used = _combined_limit_with_bin_fpr(k, fpr_for_bin, card_tag="perlifetime")
                perlifetime_r[k[2]] = r
                fh.write(f"  ctau={k[2]:g} mm:\n")
                for lb in _bin_order:
                    if lb in used:
                        fh.write(f"    {lb}: FPR = {used[lb]:g}\n")
                fh.write(f"    -> expected 95% CL limit r < {r:.4g}\n")
        print(f"Wrote per-lifetime best-FPR summary to {_pl_path}")

        # (3) Fixed-grid mode only: compare the shared combined limit across the different
        # manually-selected FPR values -- reuses the already-computed 'cells', no extra
        # combine calls (the combined-plot cells already scanned every FPR in fpr_rows).
        if not ADAPTIVE_FPR:
            _cmp_path = _LIMITS_DIR / f'fpr_comparison_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}.txt'
            with open(_cmp_path, "w") as fh:
                fh.write(f"Mass point mpi={mpi_val:g} mA={mA_val:g} "
                         f"(fixed FPR grid, shared across all lxy bins)\n")
                fh.write("Combined (all-lxy) expected limit vs FPR, per lifetime:\n")
                for k in keys:
                    fh.write(f"  ctau={k[2]:g} mm:\n")
                    for f_t, _fa, _thr in fpr_rows:
                        r = _expected_limit(cells.get((f_t, k), ()))
                        fh.write(f"    FPR={f_t:g}: expected 95% CL limit r < {r:.4g}\n")
            print(f"Wrote fixed-grid FPR comparison to {_cmp_path}")

        _best_fpr_strategy_plot(mpi_val, mA_val, keys, majority_r, perlifetime_r,
                                 fixed_cells=(cells if not ADAPTIVE_FPR else None),
                                 fpr_rows=(fpr_rows if not ADAPTIVE_FPR else None))

        _drop_key_slices()

print(f"Done. Eval plots written under {EVAL_ROOT}")
