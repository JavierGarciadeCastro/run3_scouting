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

Outputs, under BDT/<out-name>[_PNN][_cond<Vars>][_L1req]/Scenario<S>[_holdout-...]/:
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
COND_COLUMN = {"ctau": "param_ctau", "mratio": "param_mratio", "mpi": "param_mpi"}
COND_TAG    = {"ctau": "Ctau",       "mratio": "Mratio",       "mpi": "Mpi"}

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
    parser.add_argument("--scenario", choices=["A", "B1", "B2", "C"], default="A",
                        help="DQCD signal scenario to train on.")
    parser.add_argument("--sig-dir", default="tuples_DQCD_ScenarioA",
                        help="Subdirectory of the signal tuples.")
    parser.add_argument("--mass-point", default=None, metavar="mpi:mA",
                        help="Train only on this signal mass point (all its ctau values), e.g. 4:1.33.")
    parser.add_argument("--holdout", nargs="+", default=[], metavar="mpi:mA[:ctau]",
                        help="Signal point(s) excluded from training (they go to the test set only).")
    parser.add_argument("--require-l1", action=argparse.BooleanOptionalAction, default=False,
                        help="Keep only events with passL1 != 0.")

    # --- conditional / parametric inputs ---
    parser.add_argument("--conditional", action=argparse.BooleanOptionalAction, default=False,
                        help="Add the signal parameter(s) of --cond-var as input features.")
    parser.add_argument("--cond-var", choices=list(COND_COLUMN), nargs="+", default=["mratio"],
                        help="Signal parameter(s) used as inputs with --conditional.")
    parser.add_argument("--mratio-grid", type=float, nargs="+", default=[0.33, 0.10],
                        help="Nominal mA/mpi values onto which signal points are snapped (cond-var mratio).")

    # --- train/test split ---
    parser.add_argument("--do-random-splitting", action=argparse.BooleanOptionalAction, default=False,
                        help="Random stratified 70/30 split instead of the evtn-based one.")
    parser.add_argument("--train-max-events", type=int, default=15_000_000, metavar="N",
                        help="Cap on the number of training events (0 = no cap).")

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


def parse_holdout(specs):
    """'mpi:mA[:ctau]' strings -> list of (mpi, mA, ctau or None)."""
    points = []
    for spec in specs:
        values = [float(v) for v in spec.split(":")]
        points.append((values[0], values[1], values[2] if len(values) == 3 else None))
    return points


@dataclass
class Config:
    """Everything derived from the command line, in one place."""
    model_type: str
    scenario: str
    signal_dir: Path
    mass_point: tuple
    holdout: list
    require_l1: bool
    conditional: bool
    cond_vars: list          # e.g. ["mratio", "ctau"]
    cond_columns: list       # e.g. ["param_mratio", "param_ctau"]; empty if not conditional
    mratio_grid: np.ndarray
    random_splitting: bool
    train_max_events: int
    out_dir: Path
    pnn: PNNHyperparameters
    device: str


def build_config(args):
    model_tag   = "_PNN" if args.model == "pnn" else ""   # BDT keeps the workingpoint.py directory
    cond_tag    = "_cond" + "".join(COND_TAG[v] for v in args.cond_var) if args.conditional else ""
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
        scenario=args.scenario,
        signal_dir=TUPLES_BASE / args.sig_dir,
        mass_point=tuple(float(x) for x in args.mass_point.split(":")) if args.mass_point else None,
        holdout=parse_holdout(args.holdout),
        require_l1=args.require_l1,
        conditional=args.conditional,
        cond_vars=list(args.cond_var),
        cond_columns=[COND_COLUMN[v] for v in args.cond_var] if args.conditional else [],
        mratio_grid=np.array(sorted(set(args.mratio_grid)), dtype=np.float64),
        random_splitting=args.do_random_splitting,
        train_max_events=max(0, args.train_max_events),
        out_dir=HERE / (args.out_name + model_tag + cond_tag + l1_tag) / f"Scenario{args.scenario}{holdout_tag}",
        pnn=pnn,
        device=args.device,
    )


# ===========================================================================
# 4. Signal file inventory
# ===========================================================================
def snap_mratio(mpi, mA, mratio_grid):
    """mA/mpi snapped to the closest nominal value of the grid."""
    return float(mratio_grid[np.argmin(np.abs(mratio_grid - mA / mpi))])


def find_signal_files(cfg):
    """List of (path, mpi, mA, ctau), one file per signal point (merged files preferred)."""
    def to_float(s):
        return float(s.replace("p", "."))

    pattern = re.compile(rf"tuples_Signal_Scenario{cfg.scenario}_(?:(?:Par|Priv)_)?2024_"
                         rf"mpi-(\w+)_mA-(\w+)_ctau-(\w+)mm_2024(?:_\w+)?\.root")
    files_by_point = {}   # (mpi, mA, ctau) -> (path, is_merged)
    for path in sorted(glob.glob(str(cfg.signal_dir / f"tuples_Signal_Scenario{cfg.scenario}_*2024_*.root"))):
        name  = os.path.basename(path)
        match = pattern.search(name)
        if not match:
            continue
        point     = tuple(to_float(s) for s in match.groups())
        is_merged = name.endswith("_2024.root")
        previous  = files_by_point.get(point)
        if previous is None or (is_merged and not previous[1]):
            files_by_point[point] = (path, is_merged)

    signal_files = [(path, *point) for point, (path, _) in sorted(files_by_point.items())]

    if cfg.mass_point is not None:
        mpi, mA = cfg.mass_point
        signal_files = [f for f in signal_files if np.isclose(f[1], mpi) and np.isclose(f[2], mA)]
        if not signal_files:
            raise SystemExit(f"--mass-point {mpi:g}:{mA:g} matched no signal point on disk")

    if not signal_files:
        raise SystemExit(f"No signal files found in {cfg.signal_dir}")
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
    for path, mpi, mA, ctau in signal_files:
        df = read_tuple(path, cfg)
        df["param_mpi"]    = mpi
        df["param_mA"]     = mA
        df["param_ctau"]   = ctau
        df["param_mratio"] = snap_mratio(mpi, mA, cfg.mratio_grid)
        df["label"]        = 1
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_background(cfg):
    df = read_tuple(TUPLES_BASE / MINBIAS_SUBDIR / MINBIAS_FILE, cfg)
    df["label"]       = 0
    df["xsec_weight"] = MINBIAS_XSEC / MINBIAS_NGEN_AFTERFILTER
    return df


def assign_background_parameters(df_sig, df_bkg, cond_columns, seed=RANDOM_SEED):
    """Give each background event the parameter values of a random REAL signal point
    (the full combination, never a mix of values from different points)."""
    signal_points = df_sig[cond_columns].drop_duplicates().to_numpy()
    choice = np.random.default_rng(seed).integers(0, len(signal_points), size=len(df_bkg))
    df_bkg[cond_columns] = signal_points[choice]


def compute_event_weights(df):
    """Background: cross-section weight. Signal: constant, so that sum(sig) == sum(bkg)."""
    is_sig  = df["label"].to_numpy() == 1
    weights = np.ones(len(df), dtype=float)
    weights[~is_sig] = df.loc[~is_sig, "xsec_weight"].to_numpy()

    n_sig, sum_bkg = int(is_sig.sum()), float(weights[~is_sig].sum())
    if n_sig > 0 and sum_bkg > 0:
        weights[is_sig] = sum_bkg / n_sig
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

    if cfg.conditional:
        assign_background_parameters(df_sig, df_bkg, cfg.cond_columns)

    common_columns = set(df_sig.columns) & set(df_bkg.columns)
    df = pd.concat([df_sig, df_bkg], ignore_index=True)
    df["weight"] = compute_event_weights(df)
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
def holdout_mask(df, holdout):
    """True for signal events of the held-out points."""
    mask = np.zeros(len(df), dtype=bool)
    is_sig = df["label"].to_numpy() == 1
    for mpi, mA, ctau in holdout:
        sel = is_sig & np.isclose(df["param_mpi"], mpi) & np.isclose(df["param_mA"], mA)
        if ctau is not None:
            sel &= np.isclose(df["param_ctau"], ctau)
        if not sel.any():
            print(f"Holdout {mpi:g}:{mA:g}{'' if ctau is None else f':{ctau:g}'} matched no signal events.")
        mask |= sel
    return mask


def split_train_test(df, cfg):
    """70/30 split, by evtn (evtn % 10 < 7 -> train) or random stratified.
    Held-out signal points always go to test. Returns (train_rows, test_rows)."""
    is_holdout = holdout_mask(df, cfg.holdout)
    candidates = np.flatnonzero(~is_holdout)

    if cfg.random_splitting:
        labels = df["label"].to_numpy()[candidates]
        train_rows, test_rows = train_test_split(candidates, test_size=0.3, random_state=RANDOM_SEED, stratify=labels)
    else:
        is_train   = (df["evtn"].to_numpy()[candidates] % 10) < 7
        train_rows = candidates[is_train]
        test_rows  = candidates[~is_train]

    test_rows = np.concatenate([test_rows, np.flatnonzero(is_holdout)])

    # Cap the training size to bound memory
    if cfg.train_max_events and len(train_rows) > cfg.train_max_events:
        train_rows = np.random.default_rng(123).choice(train_rows, size=cfg.train_max_events, replace=False)

    return train_rows, test_rows


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


def fit_standardisation(X):
    """Per-feature mean and std (ignoring NaN); std = 1 for constant features."""
    mean = np.nanmean(X, axis=0).astype(np.float32)
    std  = np.nanstd(X, axis=0).astype(np.float32)
    std[~(std > 0)] = 1.0
    mean[~np.isfinite(mean)] = 0.0
    return mean, std


def standardise(X, mean, std):
    """(X - mean) / std, with NaN/inf set to 0 (= the mean). Returns a new float32 array."""
    Xs = (np.asarray(X, dtype=np.float32) - mean) / std
    Xs[~np.isfinite(Xs)] = 0.0
    return Xs


def class_balanced(y, w):
    """Normalise weights so that signal and background each sum to 0.5."""
    w = np.asarray(w, dtype=np.float64).copy()
    for label in (0, 1):
        w[y == label] *= 0.5 / w[y == label].sum()
    return w


class PNNClassifier:
    """Network + its input standardisation, with the same predict_proba as XGBoost."""

    def __init__(self, network, mean, std, features, hyperparameters, device="cpu"):
        self.network         = network.to(device)
        self.mean            = mean
        self.std             = std
        self.features        = features
        self.hyperparameters = hyperparameters
        self.device          = device

    def predict_logit(self, X, chunk_size=65536):
        """X must already be standardised."""
        self.network.eval()
        outputs = []
        with torch.no_grad():
            for start in range(0, len(X), chunk_size):
                batch = torch.as_tensor(X[start:start + chunk_size], device=self.device)
                outputs.append(self.network(batch).squeeze(1).cpu().numpy())
        return np.concatenate(outputs) if outputs else np.empty(0, dtype=np.float32)

    def predict_proba(self, X):
        """Raw (non-standardised) inputs -> array (n, 2) of [P(bkg), P(sig)], like XGBoost.
        For a conditional PNN, set the parameter columns to the hypothesis to test."""
        logit = self.predict_logit(standardise(np.asarray(X), self.mean, self.std))
        p_sig = 1.0 / (1.0 + np.exp(-logit))
        return np.column_stack([1.0 - p_sig, p_sig])

    def save(self, path):
        torch.save({
            "state_dict":      self.network.state_dict(),
            "mean":            self.mean,
            "std":             self.std,
            "features":        self.features,
            "hyperparameters": asdict(self.hyperparameters),
        }, path)

    @classmethod
    def load(cls, path, device="cpu"):
        saved   = torch.load(path, map_location=device, weights_only=False)
        hp      = PNNHyperparameters(**saved["hyperparameters"])
        network = build_network(len(saved["features"]), hp)
        network.load_state_dict(saved["state_dict"])
        return cls(network, saved["mean"], saved["std"], saved["features"], hp, device)


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

    Steps: split off a validation set from the TRAINING rows -> standardise inputs
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

    # 2. Standardisation (fitted on the training part only)
    mean, std = fit_standardisation(X_train)
    X_train   = standardise(X_train, mean, std)
    X_val     = standardise(X_val, mean, std)

    # 3. Class-balanced weights: sampling probabilities (train) and loss weights (validation)
    sampling_prob = class_balanced(y_train, w_train)
    w_val         = class_balanced(y_val, w_val)

    # 4. Signal parameter points (standardised), used to re-draw background parameters
    cond_idx      = [features.index(c) for c in cond_columns]
    signal_points = np.unique(X_train[np.ix_(y_train == 1, cond_idx)], axis=0) if cond_idx else None

    # 5. Model and optimiser
    classifier = PNNClassifier(build_network(len(features), hp), mean, std, features, hp, device)
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
        "tuples_subdir":      cfg.signal_dir.name,
        "lxy_bin":            lxy_label,
        "lxy_range_cm":       [LXY_BINS[i], LXY_BINS[i + 1]],
        "signal_mA":          sorted({f[2] for f in signal_files}),
        "signal_ctau_mm":     sorted({f[3] for f in signal_files}),
        "signal_mratio":      sorted({snap_mratio(f[1], f[2], cfg.mratio_grid) for f in signal_files}),
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
    print(f"Training events: {len(train_rows)}  ->  output: {cfg.out_dir}")

    # Arrays of the training rows
    feature_pos = [df.columns.get_loc(c) for c in features]
    y_train     = df["label"].to_numpy()[train_rows]
    w_train     = df["weight"].to_numpy()[train_rows]
    lxy_train   = df["lxy_bin"].to_numpy()[train_rows]

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
