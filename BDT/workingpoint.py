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
from sklearn.metrics import roc_curve, roc_auc_score
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
_parser.add_argument("--cond-var", choices=["ctau", "mratio"], default="mratio",
                     help="Which signal parameter the conditional BDT is parametrised in "
                          "(only used with --conditional). 'ctau': the lifetime [mm] "
                          "(param_ctau). 'mratio': the mass ratio mA/mpi snapped to the "
                          "nearest value of --mratio-grid (param_mratio). Output goes to "
                          "<out-name>_condCtau or <out-name>_condMratio. Default: mratio.")
_parser.add_argument("--mratio-grid", type=float, nargs="+", default=[0.33, 0.10],
                     help="Nominal mA/mpi values the signal points are snapped onto --cond-var mratio")
_parser.add_argument("--train-max-events", type=int, default=15_000_000, metavar="N",
                     help="Cap the number of TRAINING events fed to XGBoost")
_parser.add_argument("--significance", action=argparse.BooleanOptionalAction, default=True,
                     help="Run the significance evaluation after training")
_parser.add_argument("--fpr-targets", type=float, nargs="+", default=None,
                     help="Target false-positive rate(s). Default: 1e-1 1e-2 1e-3 1e-4 "
                          "(unless --fpr-scan is given; both may be combined).")
_parser.add_argument("--fpr-scan", type=float, nargs=3, default=None, metavar=("LO", "HI", "N"),
                     help="Add N log-spaced FPR targets between LO and HI (inclusive)"
                     "Merged with any explicit --fpr-targets.")
_parser.add_argument("--mass-window-rel", type=float, default=0.1,
                     help="SV1 dimuon mass window half-width as a fraction r of the nominal "
                          "signal mass mA: signal and background yields entering the "
                          "significance are counted within [mA*(1-r), mA*(1+r)]. "
                          "Default r = 0.1 (+/-10%% of mA). Set to 0 to disable (full mass range).")
_parser.add_argument("--sig-dir", default="tuples_priv_merged",
                     help="Subdirectory holding the SIGNAL tuples. ")
_parser.add_argument("--out-name", default="significance_plots_priv",
                     help="Base output directory name under BDT/")
_parser.add_argument("--require-l1", action=argparse.BooleanOptionalAction, default=False,
                     help="Keep only events with passL1 != 0 (the L1-seed decision)")
_args = _parser.parse_args()

#Explicit --fpr-targets, plus the --fpr-scan log grid
def _resolve_fpr_targets():
    t = set(_args.fpr_targets or [])
    if _args.fpr_scan is not None:
        lo, hi, n = _args.fpr_scan
        t |= {float(f"{v:.3g}") for v in np.logspace(np.log10(lo), np.log10(hi), int(n))}
    if not t:
        t = {1e-1, 1e-2, 1e-3, 1e-4}
    return sorted(t, reverse=True)

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
COND_VAR             = _args.cond_var
COND_COL             = {"ctau": "param_ctau", "mratio": "param_mratio"}[COND_VAR]
COND_TAG             = ({"ctau": "_condCtau", "mratio": "_condMratio"}[COND_VAR] if use_conditional else "")
MRATIO_GRID          = np.array(sorted(set(_args.mratio_grid)), dtype=np.float64)
TRAIN_MAX_EVENTS     = max(0, _args.train_max_events)   # 0 = no cap on training-set size
FPR_TARGETS          = _resolve_fpr_targets()
HOLDOUT              = _parse_holdout(_args.holdout)
HOLD_TAG             = ("_holdout-" + "-".join(s.replace(":", "-").replace(".", "p") for s in _args.holdout) if _args.holdout else "")
MASS_WINDOW_REL      = _args.mass_window_rel  # SV1 mass window half-width as a fraction of mA
MASS_WINDOW_ACTIVE   = MASS_WINDOW_REL > 0.0   # 0 = no window (full mass range)
REQUIRE_L1           = _args.require_l1        # keep only passL1 != 0 events (sig + bkg)
SIG_SUBDIR           = _args.sig_dir
MINBIAS_SUBDIR       = "tuples_minbias_full"
TUPLES_SUBDIR        = SIG_SUBDIR
LUMI_FB = 109.95 # 2024 Luminosity [fb^-1]
PB_TO_FB = 1.0e3
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
PLOT_LABEL = "Scouting Asymptotic Significance"
print(f"FPR targets ({len(FPR_TARGETS)}): "  + ", ".join(f"{f:g}" for f in FPR_TARGETS))

def _mass_window_halfwidth(mA):
    return MASS_WINDOW_REL * mA if MASS_WINDOW_ACTIVE else None
def _win_label_str():
    return rf'$|m_{{\mu\mu}}-m_A|<{MASS_WINDOW_REL:g}\,m_A$' if MASS_WINDOW_ACTIVE else ''

# Identify the mass ratio for each mass point
def _mratio(mpi, mA):
    return float(MRATIO_GRID[np.argmin(np.abs(MRATIO_GRID - mA / mpi))])

# Which parametrs to use for the pBDT (only when conditional is True)
def _theta_of(key):
    return key[2] if COND_VAR == "ctau" else _mratio(key[0], key[1])

# Turn FPR target into BDT score
def _wp_from_roc(fpr_arr, thr_arr, target):
    fpr_arr = np.asarray(fpr_arr, dtype=float)
    thr_arr = np.asarray(thr_arr, dtype=float)
    ok = np.flatnonzero((fpr_arr > 0.0) & (thr_arr <= 1.0))
    if len(ok) == 0:
        return None
    j = ok[int(np.argmin(np.abs(fpr_arr[ok] - target)))]
    return float(thr_arr[j]), float(fpr_arr[j])

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

def _p2f(s):
    return float(s.replace("p", "."))

def _flabel(f):
    return f"{f:g}".replace(".", "p")

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
param_grid = [(p[3], p[2], p[1]) for p in sig_file_params]  # (ctau, mA, mpi)

# Load number of generated events per signal point
with open(_HERE / "sig_ngen_cache.json") as fh:
    _ngen_by_file = json.load(fh)

# Build dict for each signal point 
SIG_NGEN = {}
for fpath, mpi_val, mA_val, ctau_val in sig_file_params:
    n = _ngen_by_file.get(os.path.basename(fpath))
    if n is None:
        print(f"no N_gen for {os.path.basename(fpath)} -- EXCLUDED from the yield table")
        continue
    SIG_NGEN[(mpi_val, mA_val, ctau_val)] = n



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
_LOAD_BRANCHES = list(dict.fromkeys(BDT_VARIABLES + ["SV1_lxy", "SV1_mass", "SV2_mass", "passL1"]))

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
lxy_bins   = [0.0, 0.2, 1.0, 2.4, 3.1, 7.0, 11.0, 16.0, 70.0]  # Match Scouting analysis
lxy_labels = ["0p0to0p2", "0p2to1p0", "1p0to2p4", "2p4to3p1", "3p1to7p0", "7p0to11p0", "11p0to16p0", "16p0to70p0"]
lxy_pretty = {lbl: rf'[{lxy_bins[i]:g}, {lxy_bins[i+1]:g}]' for i, lbl in enumerate(lxy_labels)} # For plotting legends

df_sig['lxy_bin'] = pd.cut(df_sig['SV1_lxy'], bins=lxy_bins, labels=lxy_labels, include_lowest=True)
df_bkg['lxy_bin'] = pd.cut(df_bkg['SV1_lxy'], bins=lxy_bins, labels=lxy_labels, include_lowest=True)

# Record any missing input variables
available  = set(df_sig.columns) & set(df_bkg.columns)
input_vars = [v for v in BDT_VARIABLES if v in available]
missing    = [v for v in BDT_VARIABLES if v not in available]

cond_vars = [COND_COL] if use_conditional else []

# Build the output directories
out_dir       = OUT_ROOT
_binmodel_dir = out_dir / "models" / "binned_lxy"
_fi_lxy_dir   = out_dir / "feature_importance_by_lxy"
_roc_lxy_dir  = out_dir / "ROC_by_lxy"
for _d in (out_dir, _binmodel_dir, _fi_lxy_dir, _roc_lxy_dir):
    os.makedirs(_d, exist_ok=True)

def _mp_dir(mpi_val, mA_val, lxy_label=None):
    d = out_dir / f'mpi{_flabel(mpi_val)}' / f'mA_{_flabel(mA_val)}'
    if lxy_label is not None:
        d = d / f'lxy_{lxy_label}'
    os.makedirs(d, exist_ok=True)
    return d

# ---------------------------------------------------------------------------
# Train BDT
# ---------------------------------------------------------------------------
if use_conditional:
    # Parametric BDT conditioned on one signal parameter (ctau or mA/mpi): each
    # background event gets a value drawn uniformly from the signal grid.
    _rng     = np.random.default_rng(42)
    _g_theta = np.array(sorted(set(df_sig[COND_COL].to_numpy())), dtype=np.float64)
    df_bkg[COND_COL] = _g_theta[_rng.integers(0, len(_g_theta), size=len(df_bkg))]

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
    _tr_i, _te_i = train_test_split(
        _keep_i, test_size=0.3, random_state=42, stratify=_y_all[_keep_i])
    _te_i = np.concatenate([_te_i, _hold_i])   # holdout signal -> test only
    _hp = ', '.join(f"mpi={m:g} mA={a:g}" + ('' if c is None else f" ctau={c:g}mm") for m, a, c in HOLDOUT)

else:
    _tr_i, _te_i = train_test_split(_all_i, test_size=0.3, random_state=42, stratify=_y_all)

# Cap the TRAINING rows to bound the XGBoost fit memory
if TRAIN_MAX_EVENTS and len(_tr_i) > TRAIN_MAX_EVENTS:
    _n_before = len(_tr_i)
    _tr_i = np.random.default_rng(123).choice(_tr_i, size=TRAIN_MAX_EVENTS, replace=False)

_col_pos = [df_global.columns.get_loc(c) for c in _cols] # Column positions of BDT features
y_train, y_test = _y_all[_tr_i], _y_all[_te_i] # Labels of training and testing rows
w_train, w_test = w_g[_tr_i],   w_g[_te_i] # Event weights

# Hyperparameters
_bdt_kwargs = dict(
    n_estimators=100, max_depth=3, learning_rate=0.1,
    use_label_encoder=False, eval_metric='logloss',
    tree_method='hist', n_jobs=4,
)

# ROC curve style
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
lxy_test  = _lxy_all_a[_te_i]
y_test_a  = np.asarray(y_test)
w_test_a  = np.asarray(w_test)

lxy_train = _lxy_all_a[_tr_i]
y_train_a = np.asarray(y_train)
w_train_a = np.asarray(w_train)
del _lxy_all_a

_bin_roc_trained = {}
_bin_models = {}
_bin_thr    = {}

################################
######### BDT TRAINING #########
################################
print("Training BDTs (5 minutes)")
for lxy_label in lxy_labels:
    m_te = (lxy_test == lxy_label)
    m_tr = (lxy_train == lxy_label)
    # Check that there is enough data in the bin
    if (m_tr.sum() < 20 or m_te.sum() < 10 or
            len(np.unique(y_train_a[m_tr])) < 2 or
            len(np.unique(y_test_a[m_te])) < 2):
        print(f"Skipping bin {lxy_label} " f"(train={int(m_tr.sum())}, test={int(m_te.sum())})")
        continue

    # Train BDT
    bdt_bin = XGBClassifier(**_bdt_kwargs)
    bdt_bin.fit(df_global.iloc[_tr_i[m_tr], _col_pos], y_train_a[m_tr], sample_weight=w_train_a[m_tr])

    y_score_bin = bdt_bin.predict_proba(df_global.iloc[_te_i[m_te], _col_pos])[:, 1]
    fpr_t, tpr_t, thr_t = roc_curve(y_test_a[m_te], y_score_bin, sample_weight=w_test_a[m_te])
    auc_t = roc_auc_score(y_test_a[m_te], y_score_bin, sample_weight=w_test_a[m_te])
    _bin_roc_trained[lxy_label] = (fpr_t, tpr_t, auc_t)
    _bin_models[lxy_label] = bdt_bin
    _bin_thr[lxy_label] = (fpr_t, thr_t)

    _bin_stub = f"bdt_lxy_{lxy_label}"
    _bin_model_path = _binmodel_dir / f"{_bin_stub}.json"
    bdt_bin.save_model(str(_bin_model_path))

    _bin_fi = _dump_feature_importance(
        bdt_bin,
        rf'BDT feature importance, $l_{{xy}}$ {lxy_pretty[lxy_label]} cm',
        f'feature_importance_lxy_{lxy_label}', _fi_lxy_dir)

    # Compute score thresholds per FPR target
    _bin_wp = {}
    for f_t in FPR_TARGETS:
        _wp = _wp_from_roc(fpr_t, thr_t, f_t)
        if _wp is None:
            continue
        _thr_v, _fpr_v = _wp
        _j = int(np.argmin(np.abs(fpr_t - _fpr_v)))
        _bin_wp[f"{f_t:g}"] = {
            "threshold":    _thr_v,
            "fpr_achieved": _fpr_v,
            "tpr_achieved": float(tpr_t[_j]),
        }
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
        "auc":             float(auc_t),
        "mass_window_rel": float(MASS_WINDOW_REL),
        "wp_thresholds":   _bin_wp,
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

# ROC curve plotting
for lxy_label in lxy_labels:
    roc_t = _bin_roc_trained.get(lxy_label)
    if roc_t is None:
        continue
    fpr_t, tpr_t, auc_t = roc_t
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    ax.plot(fpr_t, tpr_t, color='#d62728', linewidth=2.0, label=rf'Per-bin BDT (AUC = {auc_t:.3f})')
    _fprs = [fpr_t]
    _tprs = [tpr_t]
    _lo = _style_roc(ax, *_fprs, logy=True, tprs=_tprs)
    ax.plot([_lo, 1], [_lo, 1], 'k--', alpha=0.4, linewidth=1.0)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(rf'ROC $l_{{xy}}$ {lxy_pretty[lxy_label]} cm')
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
    ax.text(0.02, 0.97, "Preliminary", transform=ax.transAxes, fontsize=11, fontstyle="italic", fontweight="bold", va="top", ha="left")
    ax.tick_params(direction="in", top=True, right=True, which="both")
    _fout = _roc_lxy_dir / f'ROC_bylxy_{lxy_label}.png'
    fig.savefig(_fout, dpi=150, bbox_inches='tight')
    plt.close(fig)


# Score all events
print("Scoring all events (5 minutes)")
_score_cols = input_vars + cond_vars
_score_pos  = [df_global.columns.get_loc(c) for c in _score_cols]
_scores     = np.full(len(df_global), np.nan, dtype=np.float32)
_lxy_all    = df_global['lxy_bin'].to_numpy()
_CHUNK      = 2_000_000

#Do it in chunks, to reduce memory
for lxy_label, _model in _bin_models.items():
    _rows = np.flatnonzero(_lxy_all == lxy_label)
    if len(_rows) == 0:
        continue
    for _i in range(0, len(_rows), _CHUNK):
        _idx = _rows[_i:_i + _CHUNK]
        _scores[_idx] = _model.predict_proba(
            df_global.iloc[_idx, _score_pos])[:, 1]
df_global['score'] = _scores

# Pull columns we'll need repeatedly in df_global up front to avoid repeated pandas queries
_EMPTY_I  = np.empty(0, dtype=np.int64)
_SCORE_A  = df_global['score'].to_numpy()
_LABEL_A  = df_global['label'].to_numpy()
_MASS_A   = df_global['SV1_mass'].to_numpy()
_CTAU_A   = df_global['param_ctau'].to_numpy()
_THETA_A  = df_global[COND_COL].to_numpy() if use_conditional else None

_lxy_ser = df_global['lxy_bin']
if not isinstance(_lxy_ser.dtype, pd.CategoricalDtype):
    _lxy_ser = _lxy_ser.astype(pd.CategoricalDtype(categories=lxy_labels))
_LXY_CODE = _lxy_ser.cat.codes.to_numpy()
_CODE_OF  = {str(_c): _i for _i, _c in enumerate(_lxy_ser.cat.categories)}
del _lxy_ser

# Map each signal point (mpi, mA, ctau) to the row numbers of its events in df_global
_sig_i = np.flatnonzero(_LABEL_A == 1)
_MPI_S = df_global['param_mpi'].to_numpy()[_sig_i]
_MA_S  = df_global['param_mA'].to_numpy()[_sig_i]
_CT_S  = _CTAU_A[_sig_i]
_SIG_ROWS = {}
for _ct, _ma, _mp in param_grid:
    _k = (_mp, _ma, _ct)
    if _k in _SIG_ROWS:
        continue
    _SIG_ROWS[_k] = _sig_i[(_MPI_S == _k[0]) & (_MA_S == _k[1]) & (_CT_S == _k[2])]
del _MPI_S, _MA_S, _CT_S, _sig_i

_BKG_ROWS = np.flatnonzero(_LABEL_A == 0)
_BKG_MASS = _MASS_A[_BKG_ROWS]

# Group scores by lxy bin
def _split_by_bin(rows):
    out = {}
    if len(rows) == 0:
        return out
    codes  = _LXY_CODE[rows]
    scores = _SCORE_A[rows]
    for _lbl, _cd in _CODE_OF.items():
        _m = (codes == _cd)
        if _m.any():
            out[_lbl] = scores[_m]
    return out

# Per-key slices, built on first use and reused across every FPR target
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
        rows = _BKG_ROWS
        if res is not None:
            rows = rows[(_BKG_MASS >= lo) & (_BKG_MASS <= hi)]
        n_win = len(rows)
        if use_conditional and n_win:
            rows = rows[_THETA_A[rows] == theta]
        hit = (n_win, _split_by_bin(rows))
        _BKG_SLICES[ck] = hit
    return hit

def _drop_key_slices():
    _SIG_SLICES.clear()
    _BKG_SLICES.clear()

# Group ctau values by (mpi, mA) for the per-signal-point plots below.
mpi_mA_groups = {}
for ctau_val, mA_val, mpi_val in param_grid:
    mpi_mA_groups.setdefault((mpi_val, mA_val), []).append(ctau_val)

########################
###### PLOTTING ########
########################

# Discriminant (BDT score) distribution — one per (mpi, mA) PER LXY BIN
FIGSIZE     = (8.5, 6.5)
MASS_COLORS = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4", "#e377c2"]
BKG_FACE    = "#7fc7c4"
BKG_EDGE    = "#2f5f5d"
disc_bins   = np.linspace(0.0, 1.0, 51)
disc_widths = np.diff(disc_bins)
_bkg_disc_by_bin = _split_by_bin(_BKG_ROWS)
for (mpi_val, mA_val), ctau_vals in sorted(mpi_mA_groups.items()):
    _disc_sig = {c: _split_by_bin(_SIG_ROWS.get((float(mpi_val), float(mA_val), float(c)), _EMPTY_I)) for c in ctau_vals}
    for lxy_label in _bin_models:
        bkg_disc = _bkg_disc_by_bin.get(lxy_label, np.empty(0, dtype=np.float32))
        fig, ax = plt.subplots(figsize=(7.2, 5.6), constrained_layout=True)
        hb, _ = np.histogram(bkg_disc[~np.isnan(bkg_disc)], bins=disc_bins)
        if hb.sum() > 0:
            hb = hb / hb.sum()
        ax.bar(disc_bins[:-1], hb, width=disc_widths, align='edge', color=BKG_FACE, edgecolor=BKG_EDGE, linewidth=0.6, label='Background', zorder=1)

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
        ax.set_title(rf'Discriminant $m_\pi={mpi_val:g}$, $m_A={mA_val:g}$ GeV, ' rf'$l_{{xy}}$ {lxy_pretty[lxy_label]} cm')
        ax.text(0.02, 0.97, "Preliminary", transform=ax.transAxes, fontsize=11, fontstyle="italic", fontweight="bold", va="top", ha="left")
        ax.legend(loc='upper center', fontsize=9, framealpha=0.9)
        ax.tick_params(direction="in", top=True, right=True, which="both")
        _fout = (_mp_dir(mpi_val, mA_val, lxy_label) / f'Disc_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}_lxy_{lxy_label}.png')
        fig.savefig(_fout, dpi=150, bbox_inches='tight')
        plt.close(fig)
del _bkg_disc_by_bin

# Define mass windows
_MASS_WINDOW_CACHE  = {}
def _signal_mass_window(key):
    if not MASS_WINDOW_ACTIVE: #if no window is selected, scan over all the mass range
        return None
    if key in _MASS_WINDOW_CACHE:
        return _MASS_WINDOW_CACHE[key]
    mA   = key[1]
    half = _mass_window_halfwidth(mA)
    win  = (mA - half, mA + half)
    _MASS_WINDOW_CACHE[key] = (win, mA, half)
    return _MASS_WINDOW_CACHE[key]

def _counts_above_thr(scores, t):
    return int(np.count_nonzero(scores > t))

# Convert BDT score into a physical yield
def _bkg_b_at(t, key, lxy_label, theta=None):
    sc = _bkg_slices(key)[1].get(lxy_label)
    if sc is None or len(sc) == 0:
        return 0.0, 0
    denom = MINBIAS_NGEN_AFTERFILTER
    if theta is not None:
        cov = _bkg_cov.get(float(theta))
        if cov and cov > 0:
            denom = denom * cov
    n = _counts_above_thr(sc, t)
    return MINBIAS_XSEC * PB_TO_FB * LUMI_FB * (n / denom), n

def _sig_s_at(key, t, lxy_label):
    sc = _sig_slices(key)[2].get(lxy_label)
    if sc is None or len(sc) == 0:
        return None
    n_raw   = _counts_above_thr(sc, t)
    sig_eff = n_raw / SIG_NGEN[key]
    return PB_TO_FB * LUMI_FB * sig_eff, n_raw

# Conditional mode (pBDT): each background event is scored at one theta (ctau or mA/mpi)
_bkg_cov = {}
if use_conditional:
    _U = int(len(_BKG_ROWS))
    _v, _c = np.unique(_THETA_A[_BKG_ROWS], return_counts=True)
    for _vv, _cc in zip(_v, _c):
        _bkg_cov[float(_vv)] = (_cc / _U) if _U > 0 else 1.0

##############################
#### COMPUTE SIGNIFICANCE ####
##############################
if _args.significance:
    print("Computing significance and plotting (2-3 minutes)")
    fpr_targets = list(FPR_TARGETS)
    _thr_at = {}
    for f_t in fpr_targets:
        for lxy_label, (_fpr_b, _thr_b) in _bin_thr.items():
            _wp = _wp_from_roc(_fpr_b, _thr_b, f_t)
            if _wp is None:
                print(f"Lxy bin {lxy_label} has no usable cut ")
                continue
            _thr_at[(f_t, lxy_label)] = _wp
            if _wp[1] > 3.0 * f_t:
                print(f"Lxy bin {lxy_label} cannot reach FPR={f_t:g}; closest achievable is {_wp[1]:.2e}")

    fpr_rows = []
    for f_t in fpr_targets:
        _ts = [_thr_at[(f_t, l)][0] for l in _bin_thr]
        _fa = [_thr_at[(f_t, l)][1] for l in _bin_thr]
        fpr_rows.append((f_t, float(np.median(_fa)), float(np.median(_ts))))

    # Plot significance vs ctau
    _FPR_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "*", "h", "p"]
    _FPR_COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
                   "#17becf", "#e377c2", "#8c564b", "#000000", "#bcbd22",
                   "#7f7f7f", "#aec7e8"]
    _MANY_FPR = len(fpr_rows) > 4

     # Make plot labels
    def _fpr_tex(f_t):
        e = int(np.floor(np.log10(f_t)))
        m = f_t / 10.0 ** e
        if abs(m - 10.0) < 1e-9:
            e, m = e + 1, 1.0
        return rf'10^{{{e}}}' if abs(m - 1.0) < 1e-9 else rf'{m:.3g}\times 10^{{{e}}}'

    def _significance_vs_ctau_plot(cells, keys, ctaus, mpi_val, mA_val, lxy_label=None):
        if not ctaus:
            return
        order = np.argsort(np.array(ctaus, dtype=float))
        x = np.array(ctaus, dtype=float)[order]

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
        ax.text(0.0, 1.01, PLOT_LABEL, transform=ax.transAxes, ha='left', va='bottom', fontweight='bold', fontsize=13)
        ax.text(1.0, 1.01, rf'{LUMI_FB:g} fb$^{{-1}}$ (13.6 TeV, 2024)', transform=ax.transAxes, ha='right', va='bottom', fontsize=12)
        txt = [rf'Scenario {SCENARIO}', rf'$m_{{\pi_3}} = {mpi_val:g}$ GeV', rf"$m_{{A'}} = {mA_val:g}$ GeV"]
        res = _signal_mass_window(keys[0])
        if res is not None:
            (lo, hi), _c, _h = res
            txt.append(rf'Mass window: $[{lo:.2f}, {hi:.2f}]$ GeV')
        if lxy_label is not None:
            txt.append(rf'$l_{{xy}} \in {lxy_label.replace("to", "-")}$ cm')
        ax.text(0.04, 0.96, '\n'.join(txt), transform=ax.transAxes, va='top', ha='left', fontsize=11)
        leg1 = ax.legend(handles=fpr_handles, loc='upper right', framealpha=0.9,
                         fontsize=(7 if _MANY_FPR else 9), ncol=_ncol,
                         labelspacing=(0.4 if _MANY_FPR else 0.7),
                         columnspacing=0.8, handlelength=2.0,
                         title=(r'$B$ = yield (raw MC evts)' if _MANY_FPR else None),
                         title_fontsize=7)
        ax.add_artist(leg1)
        ax.tick_params(direction='in', top=True, right=True, which='both')
        fig.tight_layout()
        _fout = _mp_dir(mpi_val, mA_val, lxy_label) / f'significance_vs_ctau_mpi{_flabel(mpi_val)}_mA{_flabel(mA_val)}.png'
        fig.savefig(_fout, dpi=130)
        plt.close(fig)

    # Significance calculation
    def _compute_significance_cells(keys, lxy_only=None):
        _bins = [lxy_only] if lxy_only is not None else list(_bin_models)
        cells = {}
        for f_t, _f_ach, _thr_display in fpr_rows:
            for k in keys:
                s_tot, b_tot, saw_s = 0.0, 0.0, False
                ns_tot, nb_tot = 0, 0
                for _lb in _bins:
                    thr = _thr_at.get((f_t, _lb))
                    if thr is None:
                        continue
                    thr = thr[0]
                    s = _sig_s_at(k, thr, _lb)
                    b, n_b = _bkg_b_at(thr, k, _lb, theta=(_theta_of(k) if use_conditional else None))
                    if s is not None:
                        s, n_s = s
                        s_tot += s
                        ns_tot += n_s
                        saw_s = True
                    b_tot  += b
                    nb_tot += n_b
                s = s_tot if saw_s else None
                b = b_tot
                if s is None or s <= 0.0 or b <= 0.0:
                    cells[(f_t, k)] = (float('nan'), 0.0, s, b, ns_tot, nb_tot)
                    continue
                Z  = float(np.sqrt(2.0 * ((s + b) * np.log1p(s / b) - s)))
                p0 = float(norm.sf(Z))
                cells[(f_t, k)] = (p0, Z, s, b, ns_tot, nb_tot)
        return cells

    # Main per-mass point loop (using all the functions above)
    for (mpi_val, mA_val), ctau_vals in sorted(mpi_mA_groups.items()):
        keys = [(mpi_val, mA_val, c) for c in sorted(ctau_vals) if (mpi_val, mA_val, c) in SIG_NGEN]
        if not keys:
            continue

        ctaus = [k[2] for k in keys]
        cells = _compute_significance_cells(keys)
        _significance_vs_ctau_plot(cells, keys, ctaus, mpi_val, mA_val)

        for lxy_label in _bin_models:
            cbl = _compute_significance_cells(keys, lxy_only=lxy_label)
            _significance_vs_ctau_plot(cbl, keys, ctaus, mpi_val, mA_val, lxy_label=lxy_label)
        _drop_key_slices()
