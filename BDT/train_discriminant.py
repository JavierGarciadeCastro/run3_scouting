#!/usr/bin/env python3
"""Train the scouting signal-vs-background discriminant, one model per lxy bin.

Two kinds of discriminant are available (--model):

  bdt  XGBoost BDT. Identical to what workingpoint.py trains, written to the
       same output directory, so evaluate_workingpoint.py keeps working.

  pnn  Parametric neural network (PyTorch), inspired by ParamNN in
       BDT/input/models.py. With --conditional the signal parameters chosen
       with --cond-var are extra inputs (a true parametric NN); without it the
       network is a plain fully connected classifier.

Both use exactly the same data, selection, event weights, lxy binning and
train/test split, so their outputs can be compared directly. The test rows
are never used for training (the PNN takes its validation set for early
stopping from the TRAINING rows).

Outputs, under BDT/<out-name>[_PNN][_cond<Vars>][_L1req]/{Scenario<S>|HAHM|Scenario<S>_HAHM}[_holdout-...]/:
  models/binned_lxy/<model>_lxy_<bin>.{json,pt}   trained model of each bin
  models/binned_lxy/<model>_lxy_<bin>_manifest.json
  feature_importance_by_lxy/feature_importance_lxy_<bin>.{png,txt}
  training_curves_by_lxy/loss_lxy_<bin>.png        (PNN only)

Usage (e.g.):
    python3 BDT/train_discriminant.py --model bdt
    python3 BDT/train_discriminant.py --model pnn --conditional --cond-var mratio ctau
"""

import argparse
import copy
import glob
import json
import os
import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pandas as pd
import uproot
import xgboost
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

try:  # PyTorch is only needed for --model pnn
    import torch
except ImportError:
    torch = None


# ===========================================================================
# 1. Fixed configuration (paths, cross sections, binning, input variables)
# ===========================================================================
HERE        = Path(__file__).resolve().parent
TUPLES_BASE = Path("/ceph/cms/store/group/Run3Scouting")

MINBIAS_SUBDIR             = "tuples_minbias"
MINBIAS_FILE               = "tuples_MinBias_Fil-DoubleMuOS43_2024_2024.root"
MINBIAS_NGEN_BEFOREFILTER  = 8.31e9
MINBIAS_NGEN_AFTERFILTER   = 409318867
MINBIAS_XSEC_BEFOREFILTER  = 1.051e7  # pb
MINBIAS_XSEC               = MINBIAS_XSEC_BEFOREFILTER * (MINBIAS_NGEN_AFTERFILTER / MINBIAS_NGEN_BEFOREFILTER)  # ~= 5.18e5 pb

# Lxy binning [cm]: one independent model is trained per bin
LXY_BINS   = [0.0, 0.2, 3.1, 11.0, 70.0]
LXY_LABELS = ["0p0to0p2", "0p2to3p1", "3p1to11p0", "11p0to70p0"]
LXY_PRETTY = {label: f"[{LXY_BINS[i]:g}, {LXY_BINS[i + 1]:g}]" for i, label in enumerate(LXY_LABELS)}

# Conditional (parametric) features: --cond-var name -> dataframe column / output tag
COND_COLUMN = {"ctau": "param_ctau", "mratio": "param_mratio", "mpi": "param_mpi", "mzd": "param_mzd"}
COND_TAG    = {"ctau": "Ctau",       "mratio": "Mratio",       "mpi": "Mpi",       "mzd": "Mzd"}
MASS_PARAMS = {"dqcd": ("mpi", "mA"), "hahm": ("mzd",)}

MIN_TRAIN_EVENTS_PER_BIN = 20
RANDOM_SEED              = 42


def make_input_variables():
    """Names of the per-event input features (same list as workingpoint.py)."""
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
    variables = []
    for sv in ("SV1", "SV2"):
        variables += [f"{sv}_{s}" for s in sv_stems]
        for mu in ("mu1", "mu2"):
            variables += [f"{sv}_{mu}_{s}" for s in mu_stems]
    return variables


INPUT_VARIABLES = make_input_variables()
# Branches read from the tuples: the inputs plus what is needed for selection and binning
BRANCHES_TO_READ = list(dict.fromkeys(INPUT_VARIABLES + ["SV1_lxy", "SV1_mass", "SV2_mass", "passL1", "evtn"]))

# Plot style (same as workingpoint.py)
hep.style.use("CMS")
plt.rcParams.update({
    "font.size": 13, "axes.labelsize": 13, "axes.titlesize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 9,
    "legend.title_fontsize": 10, "axes.linewidth": 1.0, "xtick.major.size": 5,
    "ytick.major.size": 5, "xtick.minor.size": 3, "ytick.minor.size": 3,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9,
})


# ===========================================================================
# 2. Hyperparameters
# ===========================================================================
BDT_HYPERPARAMETERS = dict(
    n_estimators=100, max_depth=3, learning_rate=0.1,
    use_label_encoder=False, eval_metric="logloss",
    tree_method="hist", n_jobs=4,
)


@dataclass
class PNNHyperparameters:
    n_layers: int = 3              # hidden layers
    n_nodes: int = 64              # nodes per hidden layer
    dropout: float = 0.0
    learning_rate: float = 1e-3
    lr_decay: float = 0.9          # lr *= lr_decay whenever the training loss goes up
    batch_size: int = 1024
    epoch_size: int = 500_000      # events sampled per epoch
    max_epochs: int = 200
    min_epochs: int = 20           # never stop before this
    patience: int = 15             # stop after this many epochs without validation improvement
    tolerance: float = 1e-4        # minimum relative improvement that counts as "better"
    validation_fraction: float = 0.2
    max_validation_events: int = 1_000_000
    max_importance_events: int = 200_000


# ===========================================================================
# 3. Command line
# ===========================================================================
def parse_args():
    parser = argparse.ArgumentParser(description="Train the scouting discriminant (BDT or PNN) per lxy bin.")

    # --- which discriminant ---
    parser.add_argument("--model", choices=["bdt", "pnn"], default="bdt",
                        help="Discriminant to train: XGBoost BDT or parametric neural network.")

    # --- signal sample and selection (same meaning as in workingpoint.py) ---
    parser.add_argument("--signal", choices=["dqcd", "hahm", "both"], default="dqcd",
                        help="Signal model: 'dqcd' (mpi, mA, ctau), 'hahm' (mzd, ctau), or 'both' "
                             "(dqcd Scenario<S> + hahm, each weighted to half of the signal weight; "
                             "no --holdout/--mass-point/--sig-dir, --conditional only with --cond-var ctau). "
                             "Same as --model in evaluate_workingpoint.py.")
    parser.add_argument("--scenario", choices=["A", "B1", "B2", "C"], default="A",
                        help="DQCD signal scenario to train on.")
    parser.add_argument("--sig-dir", default=None,
                        help="Subdirectory of the signal tuples. Default: tuples_DQCD_Scenario<S> for dqcd, "
                             "tuples_hahm for hahm.")
    parser.add_argument("--mass-point", default=None, metavar="mpi:mA | mzd",
                        help="Train only on this signal mass point (all its ctau values), e.g. 4:1.33 (dqcd) or 5 (hahm).")
    parser.add_argument("--holdout", nargs="+", default=[], metavar="mpi:mA[:ctau] | mzd[:ctau]",
                        help="Signal point(s) excluded from training (they go to the test set only).")
    parser.add_argument("--require-l1", action=argparse.BooleanOptionalAction, default=False,
                        help="Keep only events with passL1 != 0.")

    # --- conditional / parametric inputs ---
    parser.add_argument("--conditional", action=argparse.BooleanOptionalAction, default=False,
                        help="Add the signal parameter(s) of --cond-var as input features.")
    parser.add_argument("--cond-var", choices=list(COND_COLUMN), nargs="+", default=None,
                        help="Signal parameter(s) used as inputs with --conditional. "
                             "Default: mratio for dqcd, mzd for hahm, ctau for both.")
    parser.add_argument("--mratio-grid", type=float, nargs="+", default=[0.33, 0.10],
                        help="Nominal mA/mpi values onto which signal points are snapped (cond-var mratio).")

    # --- train/test split ---
    parser.add_argument("--do-random-splitting", action=argparse.BooleanOptionalAction, default=False,
                        help="Random stratified 70/30 split instead of the evtn-based one.")
    parser.add_argument("--train-max-prompt-bkg", type=int, default=10_000_000, metavar="N",
                        help="Cap on the training background events of the prompt lxy bin (fixed-seed random "
                             "subset, kept events reweighted to preserve the bin's total background weight). "
                             "All signal and all other bins are used. 0 = no cap.")

    # --- output ---
    parser.add_argument("--out-name", default="significance_plots_priv",
                        help="Base output directory name under BDT/.")

    # --- PNN only ---
    pnn = parser.add_argument_group("PNN options (only used with --model pnn)")
    defaults = PNNHyperparameters()
    pnn.add_argument("--pnn-layers",     type=int,   default=defaults.n_layers)
    pnn.add_argument("--pnn-nodes",      type=int,   default=defaults.n_nodes)
    pnn.add_argument("--pnn-dropout",    type=float, default=defaults.dropout)
    pnn.add_argument("--pnn-lr",         type=float, default=defaults.learning_rate)
    pnn.add_argument("--pnn-batch-size", type=int,   default=defaults.batch_size)
    pnn.add_argument("--pnn-epoch-size", type=int,   default=defaults.epoch_size)
    pnn.add_argument("--pnn-max-epochs", type=int,   default=defaults.max_epochs)
    pnn.add_argument("--pnn-patience",   type=int,   default=defaults.patience)
    pnn.add_argument("--device", default="auto", help="PyTorch device: auto, cpu, cuda, cuda:0, ...")

    return parser.parse_args()


def parse_holdout(specs, n_mass):
    """'masses[:ctau]' strings -> list of (masses tuple, ctau or None)."""
    points = []
    for spec in specs:
        values = [float(v) for v in spec.split(":")]
        points.append((tuple(values[:n_mass]), values[n_mass] if len(values) > n_mass else None))
    return points


@dataclass
class Config:
    """Everything derived from the command line, in one place."""
    model_type: str
    signal: str
    signal_models: list
    scenario: str
    signal_dirs: dict
    mass_point: tuple
    holdout: list
    require_l1: bool
    conditional: bool
    cond_vars: list          # e.g. ["mratio", "ctau"]
    cond_columns: list       # e.g. ["param_mratio", "param_ctau"]; empty if not conditional
    mratio_grid: np.ndarray
    random_splitting: bool
    train_max_prompt_bkg: int
    out_dir: Path
    pnn: PNNHyperparameters
    device: str


def build_config(args):
    signal_models = ["dqcd", "hahm"] if args.signal == "both" else [args.signal]
    cond_vars     = list(args.cond_var or {"hahm": ["mzd"], "both": ["ctau"]}.get(args.signal, ["mratio"]))
    if args.signal == "both" and args.conditional and cond_vars != ["ctau"]:
        raise SystemExit("--signal both only supports --cond-var ctau (the mass parameters exist for one model only)")
    n_mass        = len(MASS_PARAMS[args.signal]) if args.signal != "both" else None
    default_dirs  = {"dqcd": f"tuples_DQCD_Scenario{args.scenario}", "hahm": "tuples_hahm"}
    signal_dirs   = {m: TUPLES_BASE / (args.sig_dir if args.signal != "both" and args.sig_dir else default_dirs[m])
                     for m in signal_models}
    signal_tag    = {"dqcd": f"Scenario{args.scenario}", "hahm": "HAHM", "both": f"Scenario{args.scenario}_HAHM"}[args.signal]
    model_tag   = "_PNN" if args.model == "pnn" else ""   # BDT keeps the workingpoint.py directory
    cond_tag    = "_cond" + "".join(COND_TAG[v] for v in cond_vars) if args.conditional else ""
    l1_tag      = "_L1req" if args.require_l1 else ""
    holdout_tag = ("_holdout-" + "-".join(s.replace(":", "-").replace(".", "p") for s in args.holdout)
                   if args.holdout else "")

    pnn = PNNHyperparameters(
        n_layers=args.pnn_layers, n_nodes=args.pnn_nodes, dropout=args.pnn_dropout,
        learning_rate=args.pnn_lr, batch_size=args.pnn_batch_size, epoch_size=args.pnn_epoch_size,
        max_epochs=args.pnn_max_epochs, patience=args.pnn_patience,
    )

    return Config(
        model_type=args.model,
        signal=args.signal,
        signal_models=signal_models,
        scenario=args.scenario,
        signal_dirs=signal_dirs,
        mass_point=tuple(float(x) for x in args.mass_point.split(":")) if args.mass_point else None,
        holdout=parse_holdout(args.holdout, n_mass) if args.holdout else [],
        require_l1=args.require_l1,
        conditional=args.conditional,
        cond_vars=cond_vars,
        cond_columns=[COND_COLUMN[v] for v in cond_vars] if args.conditional else [],
        mratio_grid=np.array(sorted(set(args.mratio_grid)), dtype=np.float64),
        random_splitting=args.do_random_splitting,
        train_max_prompt_bkg=max(0, args.train_max_prompt_bkg),
        out_dir=HERE / (args.out_name + model_tag + cond_tag + l1_tag) / f"{signal_tag}{holdout_tag}",
        pnn=pnn,
        device=args.device,
    )


# ===========================================================================
# 4. Signal file inventory
# ===========================================================================
def snap_mratio(mpi, mA, mratio_grid):
    """mA/mpi snapped to the closest nominal value of the grid."""
    return float(mratio_grid[np.argmin(np.abs(mratio_grid - mA / mpi))])


def signal_file_patterns(scenario):
    return {
        "dqcd": (f"tuples_Signal_Scenario{scenario}_2024_*.root",
                 re.compile(rf"tuples_Signal_Scenario{scenario}_2024_mpi-(\w+)_mA-(\w+)_ctau-(\w+)mm_2024(?:_\w+)?\.root")),
        "hahm": ("tuples_Signal_HTo2ZdTo2mu2x_MZd-*_ctau-*mm_2024*.root",
                 re.compile(r"tuples_Signal_HTo2ZdTo2mu2x_MZd-(\w+)_ctau-(\w+)mm_2024(?:_2024)?(?:_\w+)?\.root")),
    }


def find_signal_files(cfg):
    """List of (path, point, signal_model), one file per signal point (merged files preferred).
    point = masses + (ctau,): (mpi, mA, ctau) for dqcd, (mzd, ctau) for hahm."""
    def to_float(s):
        return float(s.replace("p", "."))

    patterns = signal_file_patterns(cfg.scenario)
    files_by_point = {}   # (signal_model, point) -> (path, is_merged)
    for sig_model in cfg.signal_models:
        file_glob, pattern = patterns[sig_model]
        for path in sorted(glob.glob(str(cfg.signal_dirs[sig_model] / file_glob))):
            name  = os.path.basename(path)
            match = pattern.search(name)
            if not match:
                continue
            point     = tuple(to_float(s) for s in match.groups())
            is_merged = name.endswith("_2024.root")
            previous  = files_by_point.get((sig_model, point))
            if previous is None or (is_merged and not previous[1]):
                files_by_point[(sig_model, point)] = (path, is_merged)

    signal_files = [(path, point, sig_model) for (sig_model, point), (path, _) in sorted(files_by_point.items())]

    if cfg.mass_point is not None:
        if len(cfg.mass_point) != len(MASS_PARAMS[cfg.signal]):
            raise SystemExit(f"--mass-point needs {':'.join(MASS_PARAMS[cfg.signal])} for --signal {cfg.signal}")
        signal_files = [f for f in signal_files if all(np.isclose(a, b) for a, b in zip(f[1][:-1], cfg.mass_point))]
        if not signal_files:
            raise SystemExit(f"--mass-point {':'.join(f'{v:g}' for v in cfg.mass_point)} matched no signal point on disk")

    if not signal_files:
        raise SystemExit(f"No signal files found in {', '.join(str(d) for d in cfg.signal_dirs.values())}")
    return signal_files


# ===========================================================================
# 5. Data loading
# ===========================================================================
READ_THREADS  = max(1, int(os.environ.get("WP_READ_THREADS", "4")))
READ_EXECUTOR = ThreadPoolExecutor(max_workers=READ_THREADS) if READ_THREADS > 1 else None


def read_tuple(path, cfg, step=2_000_000):
    """Read the needed branches of one ROOT file into a float32 dataframe (low peak memory)."""
    with uproot.open(path) as f:
        tree     = f["tuples"]
        branches = [b for b in BRANCHES_TO_READ if b in set(tree.keys())]
        values   = np.empty((tree.num_entries, len(branches)), dtype=np.float32)

        row = 0
        for chunk in tree.iterate(branches, library="np", step_size=step, decompression_executor=READ_EXECUTOR):
            n = len(chunk[branches[0]])
            for j, branch in enumerate(branches):
                values[row:row + n, j] = chunk[branch]
            row += n

        df = pd.DataFrame(values, columns=branches, copy=False)
        if not cfg.random_splitting:
            # evtn is only used for the train/test split; int64 because float32 is not exact above ~16.7M
            df["evtn"] = tree["evtn"].array(library="np").astype(np.int64)

    if cfg.require_l1:
        df = df[df["passL1"] > 0.5].reset_index(drop=True)
    return df


def add_dxy_over_lxy(df):
    """Muon |dxy| divided by the lifetime-weighted lxy (not stored in the tuples)."""
    for sv in ("SV1", "SV2"):
        denominator = df[f"{sv}_lxy"] * df[f"{sv}_mass"] / df[f"{sv}_ptmm"]
        denominator = np.where(denominator > 1e-9, denominator, np.float32(1e-9))
        for mu in ("mu1", "mu2"):
            df[f"{sv}_{mu}_dxy_lxy"] = np.abs(df[f"{sv}_{mu}_dxy"]) / denominator


def load_signal(cfg, signal_files):
    frames = []
    for path, point, sig_model in signal_files:
        df = read_tuple(path, cfg)
        for param, value in zip(MASS_PARAMS[sig_model], point[:-1]):
            df[f"param_{param}"] = value
        df["param_ctau"]   = point[-1]
        if sig_model == "dqcd":
            df["param_mratio"] = snap_mratio(point[0], point[1], cfg.mratio_grid)
        df["sig_model"]    = np.int8(cfg.signal_models.index(sig_model))
        df["label"]        = 1
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_background(cfg):
    df = read_tuple(TUPLES_BASE / MINBIAS_SUBDIR / MINBIAS_FILE, cfg)
    df["label"]       = 0
    df["xsec_weight"] = MINBIAS_XSEC / MINBIAS_NGEN_AFTERFILTER
    return df


def assign_background_parameters(df, background_rows, cond_columns, seed=RANDOM_SEED):
    """Give each of the given background rows the parameter values of a random REAL signal point
    (the full combination, never a mix of values from different points)."""
    signal_points = df.loc[df["label"] == 1, cond_columns].drop_duplicates().to_numpy()
    choice = np.random.default_rng(seed).integers(0, len(signal_points), size=len(background_rows))
    df.iloc[background_rows, [df.columns.get_loc(c) for c in cond_columns]] = signal_points[choice]


def compute_event_weights(df, n_signal_models):
    """Background: cross-section weight. Signal: constant per signal model, so that sum(sig) == sum(bkg)
    and each signal model gets an equal share of the signal weight."""
    is_sig  = df["label"].to_numpy() == 1
    weights = np.ones(len(df), dtype=float)
    weights[~is_sig] = df.loc[~is_sig, "xsec_weight"].to_numpy()

    sum_bkg   = float(weights[~is_sig].sum())
    sig_model = df["sig_model"].to_numpy()
    for m in range(n_signal_models):
        sel   = is_sig & (sig_model == m)
        n_sig = int(sel.sum())
        if n_sig > 0 and sum_bkg > 0:
            weights[sel] = sum_bkg / (n_signal_models * n_sig)
    return weights


def build_dataset(cfg, signal_files):
    """One dataframe with signal + background, lxy bin and event weights.
    Also returns the input variables present in BOTH signal and background."""
    print("Building dataframe (10-15 minutes)")
    df_sig = load_signal(cfg, signal_files)
    df_bkg = load_background(cfg)

    for df in (df_sig, df_bkg):
        add_dxy_over_lxy(df)
        df["lxy_bin"] = pd.cut(df["SV1_lxy"], bins=LXY_BINS, labels=LXY_LABELS, include_lowest=True)

    common_columns = set(df_sig.columns) & set(df_bkg.columns)
    df = pd.concat([df_sig, df_bkg], ignore_index=True)
    df["weight"] = compute_event_weights(df, len(cfg.signal_models))
    return df, common_columns


def feature_list(common_columns, cfg):
    """Model inputs: per-event variables present in signal AND background first, signal parameters LAST."""
    missing = [v for v in INPUT_VARIABLES if v not in common_columns]
    if missing:
        print(f"Missing input variables (skipped): {', '.join(missing)}")
    return [v for v in INPUT_VARIABLES if v in common_columns] + cfg.cond_columns


# ===========================================================================
# 6. Train/test split
# ===========================================================================
def holdout_mask(df, holdout, mass_params):
    """True for signal events of the held-out points."""
    mask = np.zeros(len(df), dtype=bool)
    is_sig = df["label"].to_numpy() == 1
    for masses, ctau in holdout:
        sel = is_sig.copy()
        for param, value in zip(mass_params, masses):
            sel &= np.isclose(df[f"param_{param}"].to_numpy(), value)
        if ctau is not None:
            sel &= np.isclose(df["param_ctau"].to_numpy(), ctau)
        if not sel.any():
            print(f"Holdout {':'.join(f'{v:g}' for v in masses)}{'' if ctau is None else f':{ctau:g}'} matched no signal events.")
        mask |= sel
    return mask


def split_train_test(df, cfg):
    """70/30 split, by evtn (evtn % 10 < 7 -> train) or random stratified.
    Held-out signal points always go to test. Returns (train_rows, test_rows)."""
    is_holdout = holdout_mask(df, cfg.holdout, MASS_PARAMS.get(cfg.signal, ()))
    candidates = np.flatnonzero(~is_holdout)

    if cfg.random_splitting:
        labels = df["label"].to_numpy()[candidates]
        train_rows, test_rows = train_test_split(candidates, test_size=0.3, random_state=RANDOM_SEED, stratify=labels)
    else:
        is_train   = (df["evtn"].to_numpy()[candidates] % 10) < 7
        train_rows = candidates[is_train]
        test_rows  = candidates[~is_train]

    test_rows = np.concatenate([test_rows, np.flatnonzero(is_holdout)])
    return train_rows, test_rows


def cap_prompt_background(df, train_rows, cfg):
    is_prompt_bkg = ((df["lxy_bin"].to_numpy()[train_rows] == LXY_LABELS[0])
                     & (df["label"].to_numpy()[train_rows] == 0))
    positions = np.flatnonzero(is_prompt_bkg)
    if not cfg.train_max_prompt_bkg or len(positions) <= cfg.train_max_prompt_bkg:
        return train_rows, 1.0
    scale = len(positions) / cfg.train_max_prompt_bkg
    drop  = np.random.default_rng(123).choice(positions, size=len(positions) - cfg.train_max_prompt_bkg, replace=False)
    print(f"Prompt-bin training background: kept {cfg.train_max_prompt_bkg:,} of {len(positions):,} "
          f"(weights x{scale:.3f}); all other training rows kept")
    return np.delete(train_rows, drop), scale


# ===========================================================================
# 7a. Discriminant: BDT
# ===========================================================================
def train_bdt(X, y, w, features):
    """Returns (model, feature_importance) with importance = XGBoost gain."""
    model = XGBClassifier(**BDT_HYPERPARAMETERS)
    model.fit(X, y, sample_weight=w)
    importance = dict(zip(features, np.asarray(model.feature_importances_, dtype=float)))
    return model, importance


# ===========================================================================
# 7b. Discriminant: parametric neural network
# ===========================================================================
def build_network(n_inputs, hp):
    """[Linear -> Dropout -> ELU] x n_layers -> Linear. Outputs a logit (sigmoid applied at prediction)."""
    layers, width = [], n_inputs
    for _ in range(hp.n_layers):
        layers += [torch.nn.Linear(width, hp.n_nodes), torch.nn.Dropout(hp.dropout), torch.nn.ELU()]
        width = hp.n_nodes
    layers.append(torch.nn.Linear(width, 1))
    return torch.nn.Sequential(*layers)


PREPROCESS_SENTINEL_MIN_ABS  = 99.0
PREPROCESS_SENTINEL_MIN_FRAC = 1e-4
PREPROCESS_LOG_TAIL_RATIO    = 20.0
PREPROCESS_CLIP              = 5.0
PREPROCESS_FIT_EVENTS        = 2_000_000


def fit_preprocessing(X, seed=RANDOM_SEED):
    """Fit the input preprocessing on (a random subset of) the training inputs:
    placeholder values (|x| >= 99 occurring exactly in >= 1e-4 of events, e.g. normChi2 = 999)
    become missing with an extra 0/1 input; long-tailed inputs (99.9% quantile of |x| > 20x its
    median) get a symmetric log, sign(x) * log(1 + |x| / s) with s the 10% quantile of |x|;
    then per-feature standardisation (ignoring missing) and clipping to +-PREPROCESS_CLIP."""
    X = np.asarray(X, dtype=np.float32)
    if len(X) > PREPROCESS_FIT_EVENTS:
        X = X[np.sort(np.random.default_rng(seed).choice(len(X), PREPROCESS_FIT_EVENTS, replace=False))]
    n_features = X.shape[1]

    sentinels = []
    for j in range(n_features):
        column = X[:, j]
        large  = column[np.isfinite(column) & (np.abs(column) >= PREPROCESS_SENTINEL_MIN_ABS)]
        values, counts = np.unique(large, return_counts=True)
        sentinels.append([float(v) for v, c in zip(values, counts) if c >= PREPROCESS_SENTINEL_MIN_FRAC * len(column)])

    prep = {"sentinels": sentinels, "log_scale": np.zeros(n_features, dtype=np.float32),
            "mean": None, "std": None, "clip": None}
    Xc = np.array(X, dtype=np.float32, copy=True)
    for j, values in enumerate(sentinels):
        if values:
            Xc[np.isin(Xc[:, j], np.asarray(values, dtype=np.float32)), j] = np.nan

    for j in range(n_features):
        magnitude = np.abs(Xc[:, j])
        magnitude = magnitude[np.isfinite(magnitude) & (magnitude > 0)]
        if len(magnitude) > 100:
            median = np.median(magnitude)
            if median > 0 and np.quantile(magnitude, 0.999) / median > PREPROCESS_LOG_TAIL_RATIO:
                prep["log_scale"][j] = np.quantile(magnitude, 0.1)

    with np.errstate(invalid="ignore"):
        for j in np.flatnonzero(prep["log_scale"] > 0):
            Xc[:, j] = np.sign(Xc[:, j]) * np.log1p(np.abs(Xc[:, j]) / prep["log_scale"][j])
        mean = np.nanmean(Xc, axis=0).astype(np.float32)
        std  = np.nanstd(Xc, axis=0).astype(np.float32)
    std[~(std > 0)] = 1.0
    mean[~np.isfinite(mean)] = 0.0
    prep.update(mean=mean, std=std, clip=PREPROCESS_CLIP)
    return prep


def preprocess(X, prep):
    """Raw inputs -> network inputs (a new float32 array): placeholders -> missing (+ 0/1 columns
    appended at the end), symmetric log, (x - mean) / std, missing/inf -> 0 (= the mean), clip."""
    X = np.asarray(X, dtype=np.float32)
    n_features = X.shape[1]
    flagged = [j for j, values in enumerate(prep["sentinels"]) if values]
    out  = np.empty((len(X), n_features + len(flagged)), dtype=np.float32)
    body = out[:, :n_features]
    body[:] = X
    for i, j in enumerate(flagged):
        hit = np.isin(body[:, j], np.asarray(prep["sentinels"][j], dtype=np.float32))
        out[:, n_features + i] = hit
        body[hit, j] = np.nan
    with np.errstate(invalid="ignore", divide="ignore"):
        for j in np.flatnonzero(prep["log_scale"] > 0):
            body[:, j] = np.sign(body[:, j]) * np.log1p(np.abs(body[:, j]) / prep["log_scale"][j])
        body -= prep["mean"]
        body /= prep["std"]
    body[~np.isfinite(body)] = 0.0
    if prep["clip"]:
        np.clip(body, -prep["clip"], prep["clip"], out=body)
    return out


def n_network_inputs(prep):
    return len(prep["mean"]) + sum(1 for values in prep["sentinels"] if values)


def class_balanced(y, w):
    """Normalise weights so that signal and background each sum to 0.5."""
    w = np.asarray(w, dtype=np.float64).copy()
    for label in (0, 1):
        w[y == label] *= 0.5 / w[y == label].sum()
    return w


class PNNClassifier:
    """Network + its input preprocessing, with the same predict_proba as XGBoost."""

    def __init__(self, network, prep, features, hyperparameters, device="cpu"):
        self.network         = network.to(device)
        self.prep            = prep
        self.features        = features
        self.hyperparameters = hyperparameters
        self.device          = device

    def predict_logit(self, X, chunk_size=65536):
        """X must already be preprocessed."""
        self.network.eval()
        outputs = []
        with torch.no_grad():
            for start in range(0, len(X), chunk_size):
                batch = torch.as_tensor(X[start:start + chunk_size], device=self.device)
                outputs.append(self.network(batch).squeeze(1).cpu().numpy())
        return np.concatenate(outputs) if outputs else np.empty(0, dtype=np.float32)

    def predict_proba(self, X):
        """Raw (non-preprocessed) inputs -> array (n, 2) of [P(bkg), P(sig)], like XGBoost.
        For a conditional PNN, set the parameter columns to the hypothesis to test."""
        logit = self.predict_logit(preprocess(X, self.prep))
        p_sig = 1.0 / (1.0 + np.exp(-logit))
        return np.column_stack([1.0 - p_sig, p_sig])

    def save(self, path):
        torch.save({
            "state_dict":      self.network.state_dict(),
            "preprocessing":   self.prep,
            "features":        self.features,
            "hyperparameters": asdict(self.hyperparameters),
        }, path)

    @classmethod
    def load(cls, path, device="cpu"):
        saved   = torch.load(path, map_location=device, weights_only=False)
        hp      = PNNHyperparameters(**saved["hyperparameters"])
        prep    = saved.get("preprocessing") or {
            "sentinels": [[] for _ in saved["features"]],
            "log_scale": np.zeros(len(saved["features"]), dtype=np.float32),
            "mean": saved["mean"], "std": saved["std"], "clip": None}
        network = build_network(n_network_inputs(prep), hp)
        network.load_state_dict(saved["state_dict"])
        return cls(network, prep, saved["features"], hp, device)


def weighted_bce(classifier, X, y, w):
    """Weighted binary cross-entropy (w already normalised to sum to 1)."""
    logit = torch.as_tensor(classifier.predict_logit(X))
    loss  = torch.nn.functional.binary_cross_entropy_with_logits(logit, torch.as_tensor(y, dtype=torch.float32), reduction="none")
    return float((loss.numpy() * w).sum())


def run_training_epoch(classifier, optimizer, X, y, sampling_prob, cond_idx, signal_points, hp, rng):
    """One epoch: sample events with probability ~ class-balanced weight, so every batch is
    ~50% signal / 50% background and the loss needs no weights. For a conditional PNN the
    background parameters are re-drawn from the real signal points in every batch."""
    classifier.network.train()
    loss_fn  = torch.nn.BCEWithLogitsLoss()
    sampled  = rng.choice(len(y), size=hp.epoch_size, replace=True, p=sampling_prob)
    losses   = []

    for start in range(0, len(sampled), hp.batch_size):
        rows    = sampled[start:start + hp.batch_size]
        X_batch = X[rows]   # fancy indexing -> copy, safe to modify
        y_batch = y[rows]

        if cond_idx:
            bkg = np.flatnonzero(y_batch == 0)
            X_batch[np.ix_(bkg, cond_idx)] = signal_points[rng.integers(0, len(signal_points), size=len(bkg))]

        X_t = torch.as_tensor(X_batch, device=classifier.device)
        y_t = torch.as_tensor(y_batch, dtype=torch.float32, device=classifier.device)

        optimizer.zero_grad()
        loss = loss_fn(classifier.network(X_t).squeeze(1), y_t)
        loss.backward()
        optimizer.step()
        losses.append(loss.item())

    return float(np.mean(losses))


def permutation_importance(classifier, X, y, w, features, rng):
    """Importance of each feature = drop in weighted ROC AUC when that column is shuffled."""
    def auc(X_):
        return roc_auc_score(y, classifier.predict_logit(X_), sample_weight=w)

    baseline   = auc(X)
    importance = {}
    for j, name in enumerate(features):
        X_perm       = X.copy()
        X_perm[:, j] = rng.permutation(X_perm[:, j])
        importance[name] = baseline - auc(X_perm)
    return importance


def train_pnn(X, y, w, features, cond_columns, hp, device):
    """Returns (classifier, feature_importance, history).

    Steps: split off a validation set from the TRAINING rows -> preprocess inputs
    (fit on train) -> train with class-balanced sampling, Adam and early stopping on
    the validation loss -> keep the best epoch -> permutation importance on validation."""
    rng = np.random.default_rng(RANDOM_SEED)
    torch.manual_seed(RANDOM_SEED)

    # 1. Training / validation split (inside the training rows)
    is_val = rng.random(len(y)) < hp.validation_fraction
    val_rows = np.flatnonzero(is_val)
    if len(val_rows) > hp.max_validation_events:
        val_rows = rng.choice(val_rows, size=hp.max_validation_events, replace=False)
    X_train, y_train, w_train = X[~is_val], y[~is_val], w[~is_val]
    X_val,   y_val,   w_val   = X[val_rows], y[val_rows], w[val_rows]

    # 2. Preprocessing (fitted on the training part only)
    prep    = fit_preprocessing(X_train)
    X_train = preprocess(X_train, prep)
    X_val   = preprocess(X_val, prep)

    # 3. Class-balanced weights: sampling probabilities (train) and loss weights (validation)
    sampling_prob = class_balanced(y_train, w_train)
    w_val         = class_balanced(y_val, w_val)

    # 4. Signal parameter points (standardised), used to re-draw background parameters
    cond_idx      = [features.index(c) for c in cond_columns]
    signal_points = np.unique(X_train[np.ix_(y_train == 1, cond_idx)], axis=0) if cond_idx else None

    # 5. Model and optimiser
    classifier = PNNClassifier(build_network(n_network_inputs(prep), hp), prep, features, hp, device)
    optimizer  = torch.optim.Adam(classifier.network.parameters(), lr=hp.learning_rate)

    # 6. Training loop with early stopping
    history     = {"train_loss": [], "val_loss": [], "lr": []}
    best_loss   = np.inf
    best_state  = None
    bad_epochs  = 0
    for epoch in range(hp.max_epochs):
        train_loss = run_training_epoch(classifier, optimizer, X_train, y_train, sampling_prob,
                                        cond_idx, signal_points, hp, rng)
        val_loss   = weighted_bce(classifier, X_val, y_val, w_val)
        lr         = optimizer.param_groups[0]["lr"]
        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)
        history["lr"].append(lr)
        print(f"    epoch {epoch + 1:3d}  train loss {train_loss:.5f}  val loss {val_loss:.5f}  lr {lr:.2e}")

        # Keep the best model
        if val_loss < best_loss * (1 - hp.tolerance):
            best_loss, best_state, bad_epochs = val_loss, copy.deepcopy(classifier.network.state_dict()), 0
        else:
            bad_epochs += 1

        # Lower the learning rate when the training loss goes up
        if len(history["train_loss"]) > 1 and history["train_loss"][-1] > history["train_loss"][-2]:
            for group in optimizer.param_groups:
                group["lr"] *= hp.lr_decay

        if epoch + 1 >= hp.min_epochs and bad_epochs >= hp.patience:
            print(f"    early stop: no validation improvement in {hp.patience} epochs")
            break

    classifier.network.load_state_dict(best_state)

    # 7. Feature importance on (a subset of) the validation set
    n_imp = min(len(y_val), hp.max_importance_events)
    importance = permutation_importance(classifier, X_val[:n_imp], y_val[:n_imp], w_val[:n_imp], features, rng)

    return classifier, importance, history


def choose_device(requested):
    if requested != "auto":
        return requested
    return "cuda" if torch.cuda.is_available() else "cpu"


# ===========================================================================
# 8. Outputs: plots and manifest
# ===========================================================================
def rank_features(importance):
    """[(name, value)] from most to least important (same tie order as workingpoint.py)."""
    names, values = list(importance), np.array(list(importance.values()), dtype=float)
    return [(names[i], values[i]) for i in np.argsort(values)[::-1]]


def plot_feature_importance(importance, title, xlabel, out_stem, top_n=30):
    """Horizontal bar plot of the top_n features (png) and the full ranking (txt)."""
    ranked = rank_features(importance)
    shown  = ranked[:top_n][::-1]   # most important on top

    fig, ax = plt.subplots(figsize=(7.5, max(4.0, 0.28 * len(shown))), constrained_layout=True)
    ax.barh(range(len(shown)), [v for _, v in shown], color="#1f77b4", edgecolor="#0f3b5f", linewidth=0.5)
    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels([name for name, _ in shown], fontsize=7)
    ax.set_xlabel(xlabel)
    ax.set_title(f"{title} (top {len(shown)})")
    ax.text(0.98, 0.02, "Preliminary", transform=ax.transAxes, fontsize=11,
            fontstyle="italic", fontweight="bold", va="bottom", ha="right")
    ax.tick_params(direction="in", top=True, right=True, which="both")
    fig.savefig(out_stem.with_suffix(".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    lines = ["rank  importance  feature"] + [f"{r:>4}  {v:>10.5f}  {name}" for r, (name, v) in enumerate(ranked)]
    out_stem.with_suffix(".txt").write_text("\n".join(lines) + "\n")


def plot_training_curves(history, title, out_path):
    epochs = np.arange(1, len(history["train_loss"]) + 1)
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.plot(epochs, history["train_loss"], label="Training")
    ax.plot(epochs, history["val_loss"], label="Validation")
    best = int(np.argmin(history["val_loss"])) + 1
    ax.axvline(best, color="gray", linestyle="--", linewidth=1, label=f"Best epoch ({best})")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Binary cross-entropy")
    ax.set_title(title)
    ax.legend()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_manifest(path, model_file, lxy_label, features, importance, cfg, signal_files, extra):
    """JSON with everything needed to reuse the model (keys as in workingpoint.py, plus model_type)."""
    i = LXY_LABELS.index(lxy_label)
    manifest = {
        "model_type":         cfg.model_type,
        "model_file":         model_file.name,
        "features":           features,
        "conditional":        cfg.conditional,
        "cond_var":           cfg.cond_vars if cfg.conditional else None,
        "cond_vars":          cfg.cond_columns,
        "require_l1":         cfg.require_l1,
        "tuples_subdir":      "+".join(d.name for d in cfg.signal_dirs.values()),
        "lxy_bin":            lxy_label,
        "lxy_range_cm":       [LXY_BINS[i], LXY_BINS[i + 1]],
        **({"signal_mA":      sorted({p[1] for _, p, m in signal_files if m == "dqcd"}),
            "signal_mratio":  sorted({snap_mratio(p[0], p[1], cfg.mratio_grid) for _, p, m in signal_files if m == "dqcd"})}
           if "dqcd" in cfg.signal_models else {}),
        **({"signal_mzd":     sorted({p[0] for _, p, m in signal_files if m == "hahm"})}
           if "hahm" in cfg.signal_models else {}),
        "signal_ctau_mm":     sorted({p[-1] for _, p, _m in signal_files}),
        "model":              cfg.signal,
        "scenario":           cfg.scenario,
        "feature_importance": {name: float(v) for name, v in rank_features(importance)},
        **extra,
    }
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)


# ===========================================================================
# 9. Main: train one model per lxy bin
# ===========================================================================
def train_one_bin(lxy_label, X, y, w, features, cfg, dirs, signal_files):
    """Train, save and document the model of one lxy bin."""
    stub  = f"{cfg.model_type}_lxy_{lxy_label}"
    title = rf"{cfg.model_type.upper()} feature importance, $l_{{xy}}$ {LXY_PRETTY[lxy_label]} cm"

    if cfg.model_type == "bdt":
        model, importance = train_bdt(X, y, w, features)
        model_file = dirs["models"] / f"{stub}.json"
        model.save_model(str(model_file))
        importance_label = "Feature importance (gain)"
        extra = {"xgboost_version": xgboost.__version__}

    else:
        model, importance, history = train_pnn(np.asarray(X, dtype=np.float32), y, w, features,
                                               cfg.cond_columns, cfg.pnn, choose_device(cfg.device))
        model_file = dirs["models"] / f"{stub}.pt"
        model.save(model_file)
        plot_training_curves(history, rf"PNN, $l_{{xy}}$ {LXY_PRETTY[lxy_label]} cm",
                             dirs["curves"] / f"loss_lxy_{lxy_label}.png")
        importance_label = "Permutation importance (AUC drop)"
        extra = {"torch_version": torch.__version__,
                 "hyperparameters": asdict(cfg.pnn),
                 "best_epoch": int(np.argmin(history["val_loss"])) + 1}

    plot_feature_importance(importance, title, importance_label,
                            dirs["importance"] / f"feature_importance_lxy_{lxy_label}")
    write_manifest(dirs["models"] / f"{stub}_manifest.json", model_file, lxy_label,
                   features, importance, cfg, signal_files, extra)


def main():
    cfg = build_config(parse_args())
    if cfg.model_type == "pnn" and torch is None:
        raise SystemExit("--model pnn needs PyTorch, which is not installed in this environment.")
    if cfg.model_type == "pnn" and not cfg.conditional:
        print("Note: --model pnn without --conditional trains a plain (non-parametric) neural network.")

    # Output directories
    dirs = {
        "models":     cfg.out_dir / "models" / "binned_lxy",
        "importance": cfg.out_dir / "feature_importance_by_lxy",
        "curves":     cfg.out_dir / "training_curves_by_lxy",
    }
    for name, d in dirs.items():
        if name != "curves" or cfg.model_type == "pnn":
            d.mkdir(parents=True, exist_ok=True)

    # Data
    signal_files = find_signal_files(cfg)
    df, common_columns = build_dataset(cfg, signal_files)
    features           = feature_list(common_columns, cfg)
    train_rows, _ = split_train_test(df, cfg)
    train_rows, prompt_bkg_scale = cap_prompt_background(df, train_rows, cfg)
    if cfg.conditional:
        assign_background_parameters(df, train_rows[df["label"].to_numpy()[train_rows] == 0], cfg.cond_columns)
    print(f"Training events: {len(train_rows)}  ->  output: {cfg.out_dir}")

    # Arrays of the training rows
    feature_pos = [df.columns.get_loc(c) for c in features]
    y_train     = df["label"].to_numpy()[train_rows]
    w_train     = df["weight"].to_numpy()[train_rows]
    lxy_train   = df["lxy_bin"].to_numpy()[train_rows]
    if prompt_bkg_scale != 1.0:
        w_train[(lxy_train == LXY_LABELS[0]) & (y_train == 0)] *= prompt_bkg_scale

    # One model per lxy bin
    print(f"Training {cfg.model_type.upper()}s")
    trained = []
    for lxy_label in LXY_LABELS:
        in_bin = lxy_train == lxy_label
        if in_bin.sum() < MIN_TRAIN_EVENTS_PER_BIN or len(np.unique(y_train[in_bin])) < 2:
            print(f"Skipping bin {lxy_label} (train events = {int(in_bin.sum())})")
            continue

        print(f"Training lxy bin {lxy_label} ({int(in_bin.sum())} events)")
        X_bin = df.iloc[train_rows[in_bin], feature_pos]
        train_one_bin(lxy_label, X_bin, y_train[in_bin], w_train[in_bin], features, cfg, dirs, signal_files)
        trained.append(lxy_label)

    if not trained:
        raise SystemExit("No lxy bin had enough events to train a model -- nothing to analyse.")
    untrained = [label for label in LXY_LABELS if label not in trained]
    if untrained:
        print(f"No model for lxy bin(s) {', '.join(untrained)} -- their events are EXCLUDED downstream.")


if __name__ == "__main__":
    main()
