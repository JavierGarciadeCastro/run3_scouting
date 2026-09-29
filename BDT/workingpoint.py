#!/usr/bin/env python3
"""BDT training and working point analysis for the scouting BDT.

Trains the global BDT, finds the BDT score threshold configurable
background rejection targets, and outputs ROCs, discriminants and
significance plots.

Usage (e.g):
    python3 BDT/workingpoint.py --fpr-scan 1e-2 1e-5 10 --out-name significance_10FPRs_nonCond
"""

import argparse
import glob
import json
import re
from concurrent.futures import ThreadPoolExecutor
import uproot
import xgboost
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mplhep as hep
import numpy as np
import pandas as pd
import os
from pathlib import Path
from sklearn.model_selection import train_test_split
from scipy.stats import norm
from xgboost import XGBClassifier

hep.style.use("CMS")
plt.rcParams.update({
    "font.size": 13, "axes.labelsize": 13, "axes.titlesize":  13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 9,
    "legend.title_fontsize": 10, "axes.linewidth": 1.0, "xtick.major.size": 5,
    "ytick.major.size": 5, "xtick.minor.size": 3, "ytick.minor.size": 3,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9,})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
_parser = argparse.ArgumentParser(description="Working point analysis for the scouting BDT.")
_parser.add_argument("--scenario", choices=["A", "B1", "B2", "C"], default="A",
                     help="DQCD signal scenario to train and evaluate on. Output goes to "
                          "significance_plots[_L1req]/Scenario<S>/.")
_parser.add_argument("--holdout", nargs="+", default=[], metavar="mpi:mA[:ctau]",
                     help="Signal point(s) to exclude from training and evaluate on the test set "
                          "only, e.g. --holdout 4:1.33 (all ctau) or --holdout 4:1.33:10.")
_parser.add_argument("--conditional", action=argparse.BooleanOptionalAction, default=False,
                     help="Train a parametric (conditional) BDT: one signal parameter (chosen "
                          "with --cond-var) is added as an input feature.")
_parser.add_argument("--cond-var", choices=["ctau", "mratio", "mpi"], default=["mratio"],
                     nargs="+",
                     help="Which signal parameter(s) the conditional BDT is parametrised in "
                          "(only used with --conditional). One or more of: 'ctau' the "
                          "lifetime [mm] (param_ctau); 'mratio' the mass ratio mA/mpi "
                          "snapped to the nearest value of --mratio-grid (param_mratio); "
                          "'mpi' the pi_3 mass [GeV] (param_mpi). All given features are "
                          "added to the BDT together. Output goes to <out-name>_cond<Tags>. "
                          "Default: mratio.")
_parser.add_argument("--mratio-grid", type=float, nargs="+", default=[0.33, 0.10],
                     help="Nominal mA/mpi values the signal points are snapped onto --cond-var mratio")
_parser.add_argument("--train-max-events", type=int, default=15_000_000, metavar="N",
                     help="Cap the number of TRAINING events fed to XGBoost")
_parser.add_argument("--sig-dir", default="tuples_priv_merged",
                     help="Subdirectory holding the SIGNAL tuples. ")
_parser.add_argument("--out-name", default="significance_plots_priv",
                     help="Base output directory name under BDT/")
_parser.add_argument("--require-l1", action=argparse.BooleanOptionalAction, default=False,
                     help="Keep only events with passL1 != 0 (the L1-seed decision)")
_parser.add_argument("--do-random-splitting", action=argparse.BooleanOptionalAction, default=False,
                     help="Use the random (stratified) train/test split instead of the evtn-based one")
_parser.add_argument("--mass-point", default=None, metavar="mpi:mA",
                     help="Train only on this single signal mass point, all its ctau values "
                          "included. E.g. --mass-point 4:1.33")
_args = _parser.parse_args()

# Convert the holdout strings into numbers
def _parse_holdout(specs):
    out = []
    for s in specs:
        parts = [float(p) for p in s.split(":")]
        mpi, mA = parts[0], parts[1]
        ctau = parts[2] if len(parts) == 3 else None
        out.append((mpi, mA, ctau))
    return out

# Set all the necessary configurations
SCENARIO             = _args.scenario
use_conditional      = _args.conditional
COND_VAR             = list(_args.cond_var)   # e.g. ["mratio"] or ["ctau", "mratio"]
_COND_COL_OF         = {"ctau": "param_ctau", "mratio": "param_mratio", "mpi": "param_mpi"}
_COND_TAG_OF         = {"ctau": "Ctau", "mratio": "Mratio", "mpi": "Mpi"}
COND_COL             = [_COND_COL_OF[v] for v in COND_VAR]   # BDT feature column(s) for the cond-var(s)
COND_TAG             = ("_cond" + "".join(_COND_TAG_OF[v] for v in COND_VAR) if use_conditional else "")
MRATIO_GRID          = np.array(sorted(set(_args.mratio_grid)), dtype=np.float64)
TRAIN_MAX_EVENTS     = max(0, _args.train_max_events)   # 0 = no cap on training-set size
HOLDOUT              = _parse_holdout(_args.holdout)
HOLD_TAG             = ("_holdout-" + "-".join(s.replace(":", "-").replace(".", "p") for s in _args.holdout) if _args.holdout else "")
REQUIRE_L1           = _args.require_l1        # keep only passL1 != 0 events (sig + bkg)
DO_RANDOM_SPLITTING  = _args.do_random_splitting  # True: random stratified split, False: evtn-based
FINAL_MASS_POINT     = (tuple(float(x) for x in _args.mass_point.split(":"))
                         if _args.mass_point else None)  # (mpi, mA) filter on the signal files
SIG_SUBDIR           = _args.sig_dir
MINBIAS_SUBDIR       = "tuples_refill_minbias"
TUPLES_SUBDIR        = SIG_SUBDIR
_HERE      = Path(__file__).resolve().parent
TUPLES_BASE    = Path("/ceph/cms/store/group/Run3Scouting")
tuples_dir     = TUPLES_BASE / SIG_SUBDIR
minbias_dir    = TUPLES_BASE / MINBIAS_SUBDIR
MINBIAS_FILE = "tuples_MinBias_Fil-DoubleMuOS43_2024_2024.root"
MINBIAS_NGEN_BEFOREFILTER = 8.31e9
MINBIAS_NGEN_AFTERFILTER = 409318867
MINBIAS_XSEC_BEFOREFILTER = 1.051e7
MINBIAS_XSEC = MINBIAS_XSEC_BEFOREFILTER * (MINBIAS_NGEN_AFTERFILTER / MINBIAS_NGEN_BEFOREFILTER)  # ~= 5.18e5 pb
OUT_ROOT = (_HERE / (_args.out_name + COND_TAG + ('_L1req' if REQUIRE_L1 else '')) / f"Scenario{SCENARIO}{HOLD_TAG}")

# Identify the mass ratio for each mass point
def _mratio(mpi, mA):
    return float(MRATIO_GRID[np.argmin(np.abs(MRATIO_GRID - mA / mpi))])

# Build boolean for the held out signal points (in case you want
# to test on a signal point that hasn't been used for training)
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
        if not sel.any():
            print(f"Holding out {mpi:g}:{mA:g}"
                  f"{'' if ctau is None else f':{ctau:g}'} matched no signal events.")
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

_SIG_GLOB = f"tuples_Signal_Scenario{SCENARIO}_*_2024_*.root"
_SIG_RE = re.compile(rf"tuples_Signal_Scenario{SCENARIO}_(?:Par|Priv)_2024_mpi-(\w+)_mA-(\w+)_ctau-(\w+)mm_2024(?:_\w+)?\.root")

# Build signal file inventory (1 file per signal point)
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

# If --mass-point was given, restrict training/scoring to that signal point (all ctau values)
if FINAL_MASS_POINT is not None:
    sig_file_params = [p for p in sig_file_params
                        if np.isclose(p[1], FINAL_MASS_POINT[0]) and np.isclose(p[2], FINAL_MASS_POINT[1])]
    if not sig_file_params:
        raise SystemExit(f"--mass-point {_args.mass_point} matched no signal point on disk")




# ----------------------------------------
# BDT variables
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
_LOAD_BRANCHES = list(dict.fromkeys(BDT_VARIABLES + ["SV1_lxy", "SV1_mass", "SV2_mass", "passL1", "evtn"]))

# Data loading helpers
_READ_THREADS  = max(1, int(os.environ.get("WP_READ_THREADS", "4")))
_READ_EXECUTOR = (ThreadPoolExecutor(max_workers=_READ_THREADS) if _READ_THREADS > 1 else None)

# Load the branches (with the least peak memory possible)
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

# (Optional: Require the events to pass L1 triggers)
def apply_l1(df, src):
    if not REQUIRE_L1:
        return df
    return df[df['passL1'] > 0.5].reset_index(drop=True)

#Build per-event weight
def compute_sample_weights(df):
    y = df['label'].values
    w = np.ones(len(df), dtype=float)
    bkg_mask = (y == 0)
    w[bkg_mask] = df.loc[bkg_mask, 'xsec_weight'].values

    #Class balance
    n_sig   = int((y == 1).sum())
    sum_bkg = float(w[bkg_mask].sum())
    if n_sig > 0 and sum_bkg > 0:
        w[y == 1] = sum_bkg / n_sig
    return w

#Add lifetime-weighted lxy (not done at tuple level)
def add_dxy_lxy(df):
    for sv in ("SV1", "SV2"):
        denom = df[f"{sv}_lxy"] * df[f"{sv}_mass"] / df[f"{sv}_ptmm"]
        denom = np.where(denom > 1e-9, denom, np.float32(1e-9))
        for mu in ("mu1", "mu2"):
            df[f"{sv}_{mu}_dxy_lxy"] = np.abs(df[f"{sv}_{mu}_dxy"]) / denom

print("Building dataframe (10-15 minutes)")
#Load all the signal into 1 dataframe
sig_frames = []
for fpath, mpi_val, mA_val, ctau_val in sig_file_params:
    df = apply_l1(read_flat(fpath), Path(fpath).name)
    df['param_ctau']   = float(ctau_val)   # conditional BDT feature (--cond-var ctau)
    df['param_mA']     = float(mA_val)
    df['param_mpi']    = float(mpi_val)
    df['param_mratio'] = _mratio(float(mpi_val), float(mA_val))   # (--cond-var mratio)
    df['label']        = 1
    sig_frames.append(df)
df_sig = pd.concat(sig_frames, ignore_index=True)
del sig_frames, df # Keep only concatenated copy in memory

# Load background dataframe
df_bkg = apply_l1(read_flat(minbias_dir / MINBIAS_FILE), MINBIAS_FILE)
df_bkg['label'] = 0
df_bkg['xsec_weight'] = MINBIAS_XSEC / MINBIAS_NGEN_AFTERFILTER

add_dxy_lxy(df_sig)
add_dxy_lxy(df_bkg)


# Lxy binning (cm)
#lxy_bins   = [0.0, 0.2, 1.0, 2.4, 3.1, 7.0, 11.0, 16.0, 70.0]  # Match Scouting analysis
#lxy_labels = ["0p0to0p2", "0p2to1p0", "1p0to2p4", "2p4to3p1", "3p1to7p0", "7p0to11p0", "11p0to16p0", "16p0to70p0"]
lxy_bins   = [0.0, 0.2, 3.1, 11.0, 70.0]  # Match Scouting analysis
lxy_labels = ["0p0to0p2", "0p2to3p1", "3p1to11p0", "11p0to70p0"]

lxy_pretty = {lbl: rf'[{lxy_bins[i]:g}, {lxy_bins[i+1]:g}]' for i, lbl in enumerate(lxy_labels)} # For plotting legends

df_sig['lxy_bin'] = pd.cut(df_sig['SV1_lxy'], bins=lxy_bins, labels=lxy_labels, include_lowest=True)
df_bkg['lxy_bin'] = pd.cut(df_bkg['SV1_lxy'], bins=lxy_bins, labels=lxy_labels, include_lowest=True)

# Record any missing input variables
available  = set(df_sig.columns) & set(df_bkg.columns)
input_vars = [v for v in BDT_VARIABLES if v in available]
missing    = [v for v in BDT_VARIABLES if v not in available]

cond_vars = list(COND_COL) if use_conditional else []

# Build the output directories
out_dir       = OUT_ROOT
_binmodel_dir = out_dir / "models" / "binned_lxy"
_fi_lxy_dir   = out_dir / "feature_importance_by_lxy"
for _d in (out_dir, _binmodel_dir, _fi_lxy_dir):
    os.makedirs(_d, exist_ok=True)


# ---------------------------------------------------------------------------
# Train BDT
# ---------------------------------------------------------------------------
if use_conditional:
    # Parametric BDT conditioned on 1+ signal parameters: each background event
    # copies the full (cond-var...) combination from a real signal point, drawn
    # uniformly -- never an unphysical mix of values from different signal points.
    _rng            = np.random.default_rng(42)
    _sig_theta_grid = df_sig[COND_COL].drop_duplicates().to_numpy()
    df_bkg[COND_COL] = _sig_theta_grid[_rng.integers(0, len(_sig_theta_grid), size=len(df_bkg))]

df_global = pd.concat([df_sig, df_bkg], ignore_index=True)
del df_sig, df_bkg

y_g   = df_global['label'] # 1: signal, 0: bkg
w_g   = compute_sample_weights(df_global) # Per-event weights
df_global['weight'] = w_g
_cols  = input_vars + cond_vars # BDT feature list
_y_all = y_g.to_numpy() # Used for stratifying into y_train/y_test

_all_i = np.arange(len(df_global)) # Row index
_hold  = _holdout_mask(df_global) # (Optional, if holdout)

# Train/Test split
if _hold.any():
    _keep_i = _all_i[~_hold]
    _hold_i = _all_i[_hold]
    if DO_RANDOM_SPLITTING:
        _tr_i, _te_i = train_test_split(
            _keep_i, test_size=0.3, random_state=42, stratify=_y_all[_keep_i])
    else:
        _tr_p, _te_p = split_by_evtn(df_global['evtn'].to_numpy()[_keep_i])   # positions within _keep_i
        _tr_i, _te_i = _keep_i[_tr_p], _keep_i[_te_p]
    _te_i = np.concatenate([_te_i, _hold_i])   # holdout signal -> test only
    _hp = ', '.join(f"mpi={m:g} mA={a:g}" + ('' if c is None else f" ctau={c:g}mm") for m, a, c in HOLDOUT)

else:
    if DO_RANDOM_SPLITTING:
        _tr_i, _te_i = train_test_split(_all_i, test_size=0.3, random_state=42, stratify=_y_all)
    else:
        _tr_i, _te_i = split_by_evtn(df_global['evtn'].to_numpy())


# Cap the TRAINING rows to bound the XGBoost fit memory
if TRAIN_MAX_EVENTS and len(_tr_i) > TRAIN_MAX_EVENTS:
    _n_before = len(_tr_i)
    _tr_i = np.random.default_rng(123).choice(_tr_i, size=TRAIN_MAX_EVENTS, replace=False)

_col_pos = [df_global.columns.get_loc(c) for c in _cols] # Column positions of BDT features
y_train = _y_all[_tr_i] # Labels of training rows
w_train = w_g[_tr_i] # Event weights

# Hyperparameters
_bdt_kwargs = dict(
    n_estimators=100, max_depth=3, learning_rate=0.1,
    use_label_encoder=False, eval_metric='logloss',
    tree_method='hist', n_jobs=4,
)

# Feature importance plot, top-N input variables
feat_names = list(input_vars + cond_vars)
def _dump_feature_importance(model, title, stub, dest_dir):
    imps  = np.asarray(model.feature_importances_, dtype=float)
    order = np.argsort(imps)[::-1]
    top_n   = min(30, len(feat_names))
    top_idx = order[:top_n][::-1]

    fig, ax = plt.subplots(figsize=(7.5, max(4.0, 0.28 * top_n)), constrained_layout=True)
    ax.barh(range(top_n), imps[top_idx], color='#1f77b4', edgecolor='#0f3b5f', linewidth=0.5)
    ax.set_yticks(range(top_n))
    ax.set_yticklabels([feat_names[i] for i in top_idx], fontsize=7)
    ax.set_xlabel('Feature importance (gain)')
    ax.set_title(f'{title} (top {top_n})')
    ax.text(0.98, 0.02, "Preliminary", transform=ax.transAxes, fontsize=11,fontstyle="italic", fontweight="bold", va="bottom", ha="right")
    ax.tick_params(direction="in", top=True, right=True, which="both")
    fig.savefig(dest_dir / f'{stub}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    _fi_lines = ["rank  importance  feature"]
    for r, i in enumerate(order):
        _fi_lines.append(f'{r:>4}  {imps[i]:>10.5f}  {feat_names[i]}')
    (dest_dir / f'{stub}.txt').write_text('\n'.join(_fi_lines) + '\n')
    return {feat_names[i]: float(imps[i]) for i in order}

# Prepare the arrays for the BDT binned in lxy
_lxy_all_a = df_global['lxy_bin'].to_numpy()
lxy_train = _lxy_all_a[_tr_i]
y_train_a = np.asarray(y_train)
w_train_a = np.asarray(w_train)
del _lxy_all_a

_bin_models = {}

################################
######### BDT TRAINING #########
################################
print("Training BDTs (5 minutes)")
for lxy_label in lxy_labels:
    m_tr = (lxy_train == lxy_label)
    # Check that there is enough data in the bin
    if m_tr.sum() < 20 or len(np.unique(y_train_a[m_tr])) < 2:
        print(f"Skipping bin {lxy_label} (train={int(m_tr.sum())})")
        continue

    # Train BDT
    bdt_bin = XGBClassifier(**_bdt_kwargs)
    bdt_bin.fit(df_global.iloc[_tr_i[m_tr], _col_pos], y_train_a[m_tr], sample_weight=w_train_a[m_tr])
    _bin_models[lxy_label] = bdt_bin

    _bin_stub = f"bdt_lxy_{lxy_label}"
    _bin_model_path = _binmodel_dir / f"{_bin_stub}.json"
    bdt_bin.save_model(str(_bin_model_path))

    _bin_fi = _dump_feature_importance(
        bdt_bin,
        rf'BDT feature importance, $l_{{xy}}$ {lxy_pretty[lxy_label]} cm',
        f'feature_importance_lxy_{lxy_label}', _fi_lxy_dir)

    # Save the BDT model
    _bin_idx = lxy_labels.index(lxy_label)
    _bin_manifest = {
        "model_file":      _bin_model_path.name,
        "features":        list(input_vars + cond_vars),
        "conditional":     bool(use_conditional),
        "cond_var":        (COND_VAR if use_conditional else None),
        "cond_vars":       list(cond_vars),
        "require_l1":      bool(REQUIRE_L1),
        "tuples_subdir":   TUPLES_SUBDIR,
        "lxy_bin":         lxy_label,
        "lxy_range_cm":    [lxy_bins[_bin_idx], lxy_bins[_bin_idx + 1]],
        "signal_mA":       sorted({float(p[2]) for p in sig_file_params}),
        "signal_ctau_mm":  sorted({float(p[3]) for p in sig_file_params}),
        "signal_mratio":   sorted({_mratio(float(p[1]), float(p[2])) for p in sig_file_params}),
        "scenario":        SCENARIO,
        "feature_importance": _bin_fi,
        "xgboost_version": xgboost.__version__,
    }
    with open(_binmodel_dir / f"{_bin_stub}_manifest.json", "w") as _mf:
        json.dump(_bin_manifest, _mf, indent=2)
    print(f"Training {lxy_label}")

if not _bin_models:
    raise SystemExit("[workingpoint] no lxy bin had enough events to train a BDT -- nothing to analyse.")
_untrained = [l for l in lxy_labels if l not in _bin_models]
if _untrained:
    print(f"No BDT for lxy bin(s) {', '.join(_untrained)} -- their events are EXCLUDED from every score, table and yield.")

