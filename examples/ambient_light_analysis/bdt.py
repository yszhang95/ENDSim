import argparse
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, classification_report
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
import ROOT
import pandas as pd
from math import *

#  feature_vecsize = 59 # DUstructure
feature_vecsize = 88 # final geometry
#  feature_vecsize = 62

# some normalization calculations
fvs = {"pen_small" : np.array([[-22, 22], [-22, 22], [-20, 90]]),
       "pen_small_equalx_short" : np.array([[-22.5, 22.5], [-22.5, 22.5], [-20, 64.5]]),
       "aframe_spacing5m" : np.array([[-23, 23], [-23, 23], [-20, 70]]),
       "aframe_spacing5m_box_wRod" : np.array([[-24, 24], [-15.5, 21.5], [-20, 66]]),
       "aframe_spacing5m_box" : np.array([[-24, 24], [-15.5, 21.5], [-20, 66]]),
       "aframe_spacing5m_DUstructure" : np.array([[-24, 24], [-9, 21.5], [-20, 66]]),
       "aframe_spacing5m_hex_DU_v2" : np.array([[-24, 24], [-9, 21.5], [-20, 66]]),
      }

def mass(V, conf="long"):
    if ("rockbed" not in conf and "equalx" not in conf and "aframe" not in conf) or ("DU" in conf):
        return (V[0][1] - V[0][0])*(V[1][1] - V[1][0])*(V[2][1] - V[2][0])*997./1.E6
    else:
        return (V[0][1] - V[0][0])*(V[2][1] - V[2][0])*((V[1][1] - (V[1][0]+10))*997. + (10)*2700.)/1.E6

frates = ROOT.TFile("/work/ENDSim/data/event_rates_lp3_depth_numu.root", "read") # I am inside docker ubuntu2404-cvmfs.
rate_perPOT_perKT = frates.Get("h_numu_lp3_me").Integral()/10.
rate_perPOT = rate_perPOT_perKT*mass(fvs["aframe_spacing5m_hex_DU_v2"], "aframe_spacing5m_hex_DU_v2")

cosmic_surface = {
                    "pen_small_equalx_short" : np.array([[-30., 30.], [-30., 100.]]),
                    "aframe_spacing5m" : np.array([[-30., 30.], [-30., 80.]]),
                    "aframe_spacing5m_box" : np.array([[-30., 30.], [-30., 80.]]),
                    "aframe_spacing5m_DUstructure" : np.array([[-30., 30.], [-30., 80.]]),
                    "aframe_spacing5m_hex_DU_v2" : np.array([[-30., 30.], [-30., 80.]])
                 }
def area(surface):
    return (surface[0][1]-surface[0][0])*(surface[1][1]-surface[1][0])
cosmic_rates = {conf: 1.5*area(cosmic_surface[conf]) for conf in cosmic_surface}
cosmic_rate = cosmic_rates["aframe_spacing5m_hex_DU_v2"]

# ---- Argument parsing -------------------------------------------------------
parser = argparse.ArgumentParser(description="Train BDT on streaming ROOT data")
parser.add_argument('--signal_root', required=True,
                    help='Path to signal ROOT file (output of pmt_bdt.cc)')
parser.add_argument('--bkg_root', required=True,
                    help='Path to background ROOT file (output of pmt_bdt.cc)')
parser.add_argument('--n_signal', type=int, required=True, default=0,
                    help='Max number of signal events to load')
parser.add_argument('--n_bkg', type=int, required=True, default=0,
                    help='Max number of background events to load')
parser.add_argument('--scale_pos_weight', type=float, default=None,
                    help='XGBoost scale_pos_weight (default: n_bkg/n_signal)')
parser.add_argument('--max_len', type=int, default=88,
                    help='Number of DOM slots per event (padding length, default: 88)')
parser.add_argument('--exclude_branches', nargs='*', default=[],
                    help='Branch names to exclude from training features (e.g. --exclude_branches dom_x dom_y)')
parser.add_argument('--batch_size', type=int, default=None,
    help='Events per mini-batch. None = load everything at once (default).')
args = parser.parse_args()

# ---- Streaming data loaders -------------------------------------------------

def get_event_row_boundaries(root_file):
    """
    Returns
    -------
    first_rows : int64 array, shape (n_unique_events,)
        Sorted row index of the first DOM row for each unique event.
    total_rows : int
    """
    rdf = ROOT.RDataFrame("doms", root_file)
    event_ids = np.array(rdf.AsNumpy(["event_id"])["event_id"])
    _, first_rows = np.unique(event_ids, return_index=True)
    return np.sort(first_rows), len(event_ids)


def load_event_batch(root_file, row_start, row_end,
                     max_len=88, exclude_branches=None):
    """
    Load a row-slice [row_start, row_end) from the 'doms' TTree, group by
    event_id, pad/flatten to a fixed-width feature matrix.

    Returns
    -------
    X         : np.float32 array, shape (n_events, n_feat_cols * max_len)
    info      : np.float32 array, shape (n_events, 8)
    feat_cols : list of str
    """
    INFO_COLS = ['rock_wgt', 'vtxX', 'vtxY', 'vtxZ', 'momX', 'momY', 'momZ', 'muE']
    if exclude_branches is None:
        exclude_branches = []

    rdf = ROOT.RDataFrame("doms", root_file).Range(int(row_start), int(row_end))
    df  = pd.DataFrame(rdf.AsNumpy())

    feat_cols = [c for c in df.columns
                 if c != 'event_id' and c not in INFO_COLS
                 and c not in exclude_branches]
    colagg = {c: list for c in feat_cols}
    colagg.update({c: 'mean' for c in INFO_COLS})
    df = df.groupby('event_id').agg(colagg).reset_index()

    for c in feat_cols:
        df[c] = df[c].apply(
            lambda l: l[:max_len] if len(l) >= max_len
                      else l + [-1.0] * (max_len - len(l))
        )
    X    = np.hstack([np.array(df[c].tolist()) for c in feat_cols]).astype(np.float32)
    info = df[INFO_COLS].to_numpy(dtype=np.float32)
    return X, info, feat_cols


def make_batch_ranges(first_rows, last_row, batch_size):
    """Return list of (row_start, row_end) covering first_rows in chunks."""
    ranges = []
    n = len(first_rows)
    for i in range(0, n, batch_size):
        r_start = first_rows[i]
        r_end   = first_rows[i + batch_size] if (i + batch_size) < n else last_row
        ranges.append((int(r_start), int(r_end)))
    return ranges


def stream_events_from_root(root_file, max_events, max_len=88, exclude_branches=None):
    """
    Thin wrapper around load_event_batch for backward-compatible full loads.

    Returns
    -------
    X         : np.float32 array, shape (n_events, n_feat_cols * max_len)
    info      : np.float32 array, shape (n_events, 8)
                columns: [rock_wgt, vtxX, vtxY, vtxZ, momX, momY, momZ, muE]
    feat_cols : list of str — the ordered feature branches actually used
    """
    first_rows, total_rows = get_event_row_boundaries(root_file)
    n = min(max_events, len(first_rows))
    if n <= 0:
        n = len(first_rows)
    row_end = int(first_rows[n]) if n < len(first_rows) else total_rows
    X, info, feat_cols = load_event_batch(root_file, 0, row_end, max_len, exclude_branches)
    print('Feature columns:', feat_cols)
    if len(X) < max_events:
        print(f"Warning: only {len(X)} events found (requested {max_events}).")
    print(f"Loaded {len(X)} events from {root_file} "
          f"({len(feat_cols)} feat cols x {max_len} slots = {len(feat_cols)*max_len} features)")
    return X[:max_events], info[:max_events], feat_cols

# ---- Normalization ROOT files -----------------------------------------------
fsignal = ROOT.TFile(args.signal_root, "read")
fbkg = ROOT.TFile(args.bkg_root, "read")
tot_signal = fsignal.Get("nevts").Integral()
# tot_bkg = fbkg.Get("nevts").Integral()
# FIXME: I use the slot for cosmic to represent ambient light.
# FIXME: tot_bkg is 10 times larger than tot_signal for training?
tot_bkg = 10 * tot_signal

# now try to normalize everything to per Day
pot = tot_signal/rate_perPOT
osc_factor = 0.73
pot_perDay = 7.92E-3
norm_signal = osc_factor*pot_perDay/pot

# livetime = tot_bkg/cosmic_rate
# spill_rate = 0.937 # Hz
# spill_length = 1.E-5 # 10us
# norm_bkg = 86400*spill_rate*spill_length/livetime

# FIXME: I use the slot for cosmic to represent ambient light
# FIXME: ambient light is 50 ms / 50us = 1000 times more than beam rate
# norm_bkg = 1000*norm_signal/osc_factor
# norm_bkg = 1000*norm_signal/osc_factor
norm_bkg = 2000*norm_signal/osc_factor

print("Normalization Factors : ", norm_signal, norm_bkg)

outfolder = "end/efficiency/domPE/bdt_aframe_spacing5m_hex_DU_v2/3dom_3pe_50/"

# Load dataset by streaming from ROOT files
if args.batch_size is None:
    X_signal, df_signal, feat_cols = stream_events_from_root(args.signal_root, args.n_signal, args.max_len, args.exclude_branches)
    X_bkg,    df_bkg,    _         = stream_events_from_root(args.bkg_root,    args.n_bkg,    args.max_len, args.exclude_branches)
    X_signal[np.isinf(X_signal)] = -5
    X_bkg[np.isinf(X_bkg)] = -5

    X_info = np.vstack([df_signal, df_bkg])
    signal_eff_corr = tot_signal/(df_signal[:,0].sum())
    bkg_eff_corr = tot_bkg/(df_bkg[:,0].sum())
    print("No Activity Efficiency Correction : ", signal_eff_corr, bkg_eff_corr)
    del df_signal, df_bkg

    # Build label vectors
    y_signal = np.ones(len(X_signal),  dtype=int)
    y_bkg = np.zeros(len(X_bkg), dtype=int)
    print("Signal vs Background Statistics : ", len(y_signal), len(y_bkg))

    # Stack together
    X = np.vstack([X_signal, X_bkg])
    y = np.concatenate([y_signal, y_bkg])
    del X_signal, X_bkg
else:
    sig_first_rows, sig_total_rows = get_event_row_boundaries(args.signal_root)
    bkg_first_rows, bkg_total_rows = get_event_row_boundaries(args.bkg_root)

    # Limit to requested event count; compute the exclusive end row carefully
    if args.n_signal == 0:
        args.n_signal = len(sig_first_rows)
    if args.n_bkg == 0:
        args.n_bkg = len(bkg_first_rows)
    n_sig = min(args.n_signal, len(sig_first_rows))
    n_bkg = min(args.n_bkg,    len(bkg_first_rows))
    sig_last_row = (int(sig_first_rows[n_sig]) if n_sig < len(sig_first_rows) else sig_total_rows)
    bkg_last_row = (int(bkg_first_rows[n_bkg]) if n_bkg < len(bkg_first_rows) else bkg_total_rows)
    sig_first_rows = sig_first_rows[:n_sig]
    bkg_first_rows = bkg_first_rows[:n_bkg]

    # Validation set (first 20% of each class, loaded once)
    n_val_sig = max(1, int(n_sig * 0.2))
    n_val_bkg = max(1, int(n_bkg * 0.2))

    val_row_end_sig = (int(sig_first_rows[n_val_sig])
                       if n_val_sig < len(sig_first_rows) else sig_last_row)
    val_row_end_bkg = (int(bkg_first_rows[n_val_bkg])
                       if n_val_bkg < len(bkg_first_rows) else bkg_last_row)

    X_val_sig, info_val_sig, feat_cols = load_event_batch(
        args.signal_root, 0, val_row_end_sig, args.max_len, args.exclude_branches)
    X_val_bkg, info_val_bkg, _         = load_event_batch(
        args.bkg_root,    0, val_row_end_bkg, args.max_len, args.exclude_branches)
    X_val_sig[np.isinf(X_val_sig)] = -5
    X_val_bkg[np.isinf(X_val_bkg)] = -5

    X_val  = np.vstack([X_val_sig, X_val_bkg])
    y_val  = np.concatenate([np.ones(len(X_val_sig), dtype=int),
                              np.zeros(len(X_val_bkg), dtype=int)])
    signal_eff_corr = tot_signal / info_val_sig[:,0].sum()
    bkg_eff_corr    = tot_bkg    / info_val_bkg[:,0].sum()
    del X_val_sig, X_val_bkg, info_val_sig, info_val_bkg

    # Split: 20% val | 60% train | 20% test
    n_test_sig = max(1, int(n_sig * 0.2))
    n_test_bkg = max(1, int(n_bkg * 0.2))
    n_train_end_sig = n_sig - n_test_sig
    n_train_end_bkg = n_bkg - n_test_bkg
    train_row_end_sig = (int(sig_first_rows[n_train_end_sig])
                         if n_train_end_sig < len(sig_first_rows) else sig_last_row)
    train_row_end_bkg = (int(bkg_first_rows[n_train_end_bkg])
                         if n_train_end_bkg < len(bkg_first_rows) else bkg_last_row)

    # Training batch ranges (middle 60%)
    sig_train_ranges = make_batch_ranges(
        sig_first_rows[n_val_sig:n_train_end_sig], train_row_end_sig, args.batch_size)
    bkg_train_ranges = make_batch_ranges(
        bkg_first_rows[n_val_bkg:n_train_end_bkg], train_row_end_bkg, args.batch_size)

    # Test batch ranges (last 20%, held out from training)
    sig_test_ranges = make_batch_ranges(
        sig_first_rows[n_train_end_sig:], sig_last_row, args.batch_size)
    bkg_test_ranges = make_batch_ranges(
        bkg_first_rows[n_train_end_bkg:], bkg_last_row, args.batch_size)

    # scale_pos_weight proxy from event counts
    if args.scale_pos_weight is None:
        scale_pos_weight = n_bkg / max(n_sig, 1)
    else:
        scale_pos_weight = args.scale_pos_weight

    print(f"Batch mode: {n_sig} signal / {n_bkg} bkg events, "
          f"{len(sig_train_ranges)} train batches, {len(sig_test_ranges)} test batches")

signal_evtrate_factor = 1.
bkg_evtrate_factor = 1.

def apply_trigger(X, trigger_branch, feat_cols, max_len):
    #  # older tree
    #  ndoms_pecut = np.sum(X[:, 3*feature_vecsize:4*feature_vecsize] >= 3, axis=1)
    #  # newer tree, 100ns for pe
    #  ndoms_pecut = np.sum(X[:, 6*feature_vecsize:7*feature_vecsize] >= 3, axis=1)
    # newer tree, 50ns for pe — use dynamic index in case branches were excluded
    if trigger_branch not in feat_cols:
        raise ValueError(f"Trigger branch '{trigger_branch}' was excluded; cannot apply trigger.")
    idx = feat_cols.index(trigger_branch)
    ndoms_pecut = np.sum(X[:, idx*max_len:(idx+1)*max_len] >= 3, axis=1)
    return np.where(ndoms_pecut >= 3)

if args.batch_size is None:
    # Split: 60% train, 20% validation, 20% test
    X_train, X_temp, y_train, y_temp, info_train, info_temp = train_test_split(
        X, y, X_info, test_size=0.4, stratify=y, random_state=42
    )
    X_val, X_test, y_val, y_test, info_val, info_test = train_test_split(
        X_temp, y_temp, info_temp, test_size=0.5, stratify=y_temp, random_state=43
    )
    del X, X_info  # full arrays no longer needed

    # Create DMatrix objects
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dval   = xgb.DMatrix(X_val,   label=y_val)
    dtest  = xgb.DMatrix(X_test,  label=y_test)
    trigger_cut_nb = apply_trigger(X_test, 'npe_50', feat_cols, args.max_len)
    del X_train, X_val, X_test

    # Compute scale_pos_weight
    if args.scale_pos_weight is not None:
        scale_pos_weight = args.scale_pos_weight
    else:
        scale_pos_weight = len(y_bkg) / max(len(y_signal), 1)
else:
    # Validation DMatrix already built from the validation slice
    dval = xgb.DMatrix(X_val, label=y_val)
    del X_val  # free memory; dval is already constructed

print("scale_pos_weight : ", scale_pos_weight)

# XGBoost parameters
params = {
    'objective': 'binary:logistic',
    'eval_metric': ['auc', 'logloss'],
    'tree_method': 'hist',
    'max_depth': 7,
    'eta': 0.5,
    'subsample': 0.8,
    'colsample_bytree': 0.8,
    'seed': 30,
    'scale_pos_weight': scale_pos_weight
}

# Train
start = time.time()
if args.batch_size is None:
    # Train with early stopping on validation set
    evals = [(dtrain, 'train'), (dval, 'val')]
    bst = xgb.train(
        params,
        dtrain,
        num_boost_round=500,
        evals=evals,
        early_stopping_rounds=20,
        verbose_eval=10
    )
    print("Time Taken for Training (s) : ", time.time()-start)
else:
    n_train_batches = max(len(sig_train_ranges), len(bkg_train_ranges))
    rounds_per_batch = max(1, 500 // n_train_batches)
    print(f"Batched training: {n_train_batches} batches x {rounds_per_batch} rounds")

    bst = None
    for i in range(n_train_batches):
        sig_r = sig_train_ranges[i % len(sig_train_ranges)]
        bkg_r = bkg_train_ranges[i % len(bkg_train_ranges)]

        X_s, _, _ = load_event_batch(args.signal_root, *sig_r, args.max_len, args.exclude_branches)
        X_b, _, _ = load_event_batch(args.bkg_root,    *bkg_r, args.max_len, args.exclude_branches)
        X_s[np.isinf(X_s)] = -5
        X_b[np.isinf(X_b)] = -5

        X_batch = np.vstack([X_s, X_b])
        y_batch = np.concatenate([np.ones(len(X_s), dtype=int),
                                   np.zeros(len(X_b), dtype=int)])
        del X_s, X_b

        dtrain_batch = xgb.DMatrix(X_batch, label=y_batch)
        del X_batch

        bst = xgb.train(
            params,
            dtrain_batch,
            num_boost_round=rounds_per_batch,
            evals=[(dtrain_batch, f'train_b{i}'), (dval, 'val')],
            xgb_model=bst,
            verbose_eval=False
        )
        print(f"Batch {i+1}/{n_train_batches} done.")

    print("Time Taken for Training (s):", time.time() - start)

# Save model for later inference
bst.save_model('xgb_model_aframe_spacing5m_hex_DU_v2.json')

# Evaluate on the test set
if args.batch_size is None:
    y_pred_prob = bst.predict(dtest)
else:
    # Predict only on the held-out test set (last 20%)
    trigger_col_idx = feat_cols.index('npe_50')

    all_pred, all_y, all_info, all_trig = [], [], [], []

    def _predict_loop(root_file, ranges, label_val):
        for r in ranges:
            X_b, info_b, _ = load_event_batch(
                root_file, *r, args.max_len, args.exclude_branches)
            X_b[np.isinf(X_b)] = -5
            # Trigger mask computed while X_b is still in memory
            ndoms = np.sum(
                X_b[:, trigger_col_idx*args.max_len:(trigger_col_idx+1)*args.max_len] >= 3,
                axis=1)
            all_trig.append(ndoms >= 3)
            all_pred.append(bst.predict(xgb.DMatrix(X_b)))
            all_y.append(np.full(len(X_b), label_val, dtype=int))
            all_info.append(info_b)
            del X_b, info_b

    _predict_loop(args.signal_root, sig_test_ranges, 1)
    _predict_loop(args.bkg_root,    bkg_test_ranges, 0)

    y_pred_prob  = np.concatenate(all_pred)
    y_test       = np.concatenate(all_y)
    info_test    = np.vstack(all_info)
    trigger_mask = np.concatenate(all_trig)
    del all_pred, all_y, all_info, all_trig

    # Recompute efficiency corrections with all weights
    signal_eff_corr = tot_signal / (info_test[:,0][y_test == 1]).sum()
    bkg_eff_corr    = tot_bkg    / (info_test[:,0][y_test == 0]).sum()

y_pred = (y_pred_prob >= 0.8).astype(int)

# Print metrics
auc = roc_auc_score(y_test, y_pred_prob)
print(f"Test ROC-AUC: {auc:.4f}")
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print("Classification Report:")
print(classification_report(y_test, y_pred))

# Plot and save ROC curve
fpr, tpr, _ = roc_curve(y_test, y_pred_prob)
np.save('aframe_spacing5m_hex_DU_v2_fpr.npy', fpr)
np.save('aframe_spacing5m_hex_DU_v2_tpr.npy', tpr)
plt.figure()
plt.plot(fpr, tpr)
plt.plot([0, 1], [0, 1], '--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.grid(True)
plt.savefig(outfolder+'roc_curve.pdf')

def angle(t_info):
    mom_magnitude = np.linalg.norm(t_info[:,4:7], axis=1)
    cos_angle = np.clip(t_info[:,6]/mom_magnitude, -1.0, 1.0)
    return np.arccos(cos_angle)*180.0/pi

info_angle = angle(info_test)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
hist1 = ax1.hist2d(y_pred_prob[y_test==1], info_angle[y_test==1], bins=30, cmap='viridis', norm=colors.LogNorm())
ax1.set_title('Beam')
hist2 = ax2.hist2d(y_pred_prob[y_test==0], info_angle[y_test==0], bins=30, cmap='viridis', norm=colors.LogNorm())
ax2.set_title('Cosmics')
ax1.set_xlabel('xgboost Score')
ax1.set_ylabel('Muon Angle (Degrees)')
ax2.set_xlabel('xgboost Score')
ax2.set_ylabel('Muon Angle (Degrees)')
plt.colorbar(hist1[3], ax=ax1)
plt.colorbar(hist2[3], ax=ax2)
plt.savefig(outfolder+'bdt_anglecorr.pdf')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
hist1 = ax1.hist2d(y_pred_prob[y_test==1], info_test[:,7][y_test==1], bins=30, cmap='viridis', norm=colors.LogNorm())
ax1.set_title('Beam')
hist2 = ax2.hist2d(y_pred_prob[y_test==0], info_test[:,7][y_test==0], bins=30, cmap='viridis', norm=colors.LogNorm())
ax2.set_title('Cosmics')
ax1.set_xlabel('xgboost Score')
ax1.set_ylabel('Muon KE (GeV)')
ax2.set_xlabel('xgboost Score')
ax2.set_ylabel('Muon KE (GeV)')
plt.colorbar(hist1[3], ax=ax1)
plt.colorbar(hist2[3], ax=ax2)
plt.savefig(outfolder+'bdt_energycorr.pdf')

if args.batch_size is None:
    trigger_cut = trigger_cut_nb
else:
    trigger_cut = np.where(trigger_mask)
y_pred_prob = y_pred_prob[trigger_cut]
info_test = info_test[trigger_cut]
y_test = y_test[trigger_cut]
signal_eff_corr = tot_signal/(info_test[:, 0][y_test == 1]).sum()
bkg_eff_corr = tot_bkg/(info_test[:, 0][y_test == 0]).sum()
print("Trigger Efficiency Correction : ", signal_eff_corr, bkg_eff_corr)

signal_probs = y_pred_prob[y_test == 1]
bg_probs = y_pred_prob[y_test == 0]
signal_wgts = signal_evtrate_factor*norm_signal*info_test[:, 0][y_test == 1]
bkg_wgts = bkg_evtrate_factor*norm_bkg*info_test[:, 0][y_test == 0]
plt.figure()
plt.hist(signal_probs, bins=50, alpha=0.5, label='Beam', weights=signal_wgts)
plt.hist(signal_probs, bins=50, histtype='step', cumulative=-1, weights=signal_wgts, label='Cumulative Signal')
plt.hist(bg_probs, bins=50, alpha=0.5, label='Cosmics', weights=bkg_wgts)
plt.xlabel('xgboost Score')
plt.ylabel('Counts per Day')
plt.yscale('log')
plt.legend()
plt.grid(True)
plt.savefig(outfolder+'prob_hist.pdf')
print('min/max score for signal: ', signal_probs.min(), signal_probs.max())
print('min/max score for bkg: ', bg_probs.min(), bg_probs.max())
np.savez("prob_array.npz", bg_probs=bg_probs, signal_probs=signal_probs)

#  # get an idea of s/n over fiducial volume
#  y_cut = 4.0 #m
#  x_cuts = np.linspace(1, 10, 10)
#
#  fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
#  for x_cut in x_cuts:
#      inFV = (np.abs(info_test[:, 1]) <= 1000.0*x_cut) & (info_test[:, 2] <= 1000.0*y_cut)
#      signal_inFV = y_pred_prob[inFV & (y_test == 1)]
#      bkg_inFV = y_pred_prob[inFV & (y_test == 0)]
#      signal_inFV_wgts = signal_evtrate_factor*norm_signal*info_test[:, 0][inFV & (y_test == 1)]
#      bkg_inFV_wgts = bkg_evtrate_factor*norm_bkg*info_test[:, 0][inFV & (y_test == 0)]
#
#      # bin them
#      signal_hist, bin_edges = np.histogram(signal_inFV, bins=5, weights=signal_inFV_wgts)
#      bkg_hist, _ = np.histogram(bkg_inFV, bins=5, weights=bkg_inFV_wgts)
#      ratio = np.divide(signal_hist, signal_hist+bkg_hist, where=(signal_hist+bkg_hist != 0))
#      ax1.hist(bin_edges[:-1], bin_edges, weights=ratio, label='Canopy Width : %d m'%(2*x_cut), histtype='step')
#      ax1.set_xlabel('xgboost Score')
#      ax1.set_ylabel('S / (S + B)')
#      ax1.legend(loc='upper left')
#
#      ax2.hist(signal_inFV, bins=50, weights=signal_inFV_wgts, cumulative=-1, label='Canopy Width: %d m'%(2*x_cut), histtype='step')
#      ax2.axhline(y=2, color='red', linestyle='--', label='Program Metric')
#      ax2.set_xlabel('xgboost Score')
#      ax2.set_ylabel('Cumulative Signal Counts per Day')
#  plt.title('Canopy Height at y=4m Always (Only Width is Varied)')
#  plt.savefig(outfolder+'prob_canopy.pdf')
#
#  x_cut = 4.0 #m
#  y_floors = np.linspace(-3, -13, 11)
#  fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
#  for y_floor_cut in y_floors:
#      inFV = (info_test[:, 2] >= 1000.0*y_floor_cut) #& (info_test[:, 2] <= 1000.0*y_cut) & (np.abs(info_test[:, 1]) <= 1000.0*x_cut) # assume 6m canopy width
#      signal_inFV = y_pred_prob[inFV & (y_test == 1)]
#      bkg_inFV = y_pred_prob[inFV & (y_test == 0)]
#      signal_inFV_wgts = signal_evtrate_factor*norm_signal*info_test[:, 0][inFV & (y_test == 1)]
#      bkg_inFV_wgts = bkg_evtrate_factor*norm_bkg*info_test[:, 0][inFV & (y_test == 0)]
#
#      # bin them
#      signal_hist, bin_edges = np.histogram(signal_inFV, bins=5, weights=signal_inFV_wgts)
#      bkg_hist, _ = np.histogram(bkg_inFV, bins=5, weights=bkg_inFV_wgts)
#      ratio = np.divide(signal_hist, signal_hist+bkg_hist, where=(signal_hist+bkg_hist != 0))
#      ax1.hist(bin_edges[:-1], bin_edges, weights=ratio, label='Floor Height : %d m'%(-1*y_floor_cut - 3), histtype='step')
#      ax1.set_xlabel('xgboost Score')
#      ax1.set_ylabel('S / (S + B)')
#      #  ax1.legend(loc='upper left')
#
#      ax2.hist(signal_inFV, bins=50, weights=signal_inFV_wgts, cumulative=-1, label='Floor Height: %d m'%(-1*y_floor_cut - 3), histtype='step')
#      ax2.set_xlabel('xgboost Score')
#      ax2.set_ylabel('Cumulative Signal Counts per Day')
#      ax2.legend(loc='lower center')
#      ax2.axhline(y=2, color='red', linestyle='--')
#  #  plt.title('Canopy Width at x=(-4, 4)m and Height at y=4m')
#  plt.title('No Canopy')
#  plt.savefig(outfolder+'prob_floor_nocanopy.pdf')

thresholds = np.linspace(0, 1, 101)
eff = []
pur = []
bkg_eff = []
for thr in thresholds:
    preds = (y_pred_prob >= thr).astype(int)
    # Apply rock wgt correctly
    tp = (info_test[:,0][((preds == 1) & (y_test == 1))].sum())*norm_signal
    fp = (info_test[:,0][((preds == 1) & (y_test == 0))].sum())*norm_bkg
    fn = (info_test[:,0][((preds == 0) & (y_test == 1))].sum())*norm_signal
    tn = (info_test[:,0][((preds == 0) & (y_test == 0))].sum())*norm_bkg
    eff.append(tp / ((tp + fn)*signal_eff_corr) if (tp + fn) > 0 else 0)
    pur.append(tp / (tp + fp) if (tp + fp) > 0 else 0)
    bkg_eff.append(fp / ((fp + tn)*bkg_eff_corr) if (fp + tn) > 0 else 0)

eff = np.array(eff)
pur = np.array(pur)
bkg_eff = np.array(bkg_eff)
ep = eff * pur

plt.figure()
plt.plot(thresholds, 3*eff, label='Signal Efficiency (x3)')
plt.plot(thresholds, pur, label='Purity')
plt.plot(thresholds, 10*bkg_eff, label='Bkg. Efficiency (x10)')
plt.plot(thresholds, 3*ep,  label='Signal Eff*Pur (x3)')
plt.xlabel('xgboost Score')
plt.ylabel('Metric')
plt.legend()
plt.grid(True)
plt.savefig(outfolder+'fom_opt.pdf')

import re
importance = bst.get_score(importance_type='gain')  # or 'weight', 'cover'
keys = list(importance.keys())
values = list(importance.values())
indices = np.argsort(values)[::-1]
#  features = ['dom_x', 'dom_y', 'dom_z', 'npe', 'npmts', 'pe_min', 'pe_rms', 'pe_spread', 't_mean', 't_min', 't_rms', 't_spread']
#  keys_names = []
#  for i in indices:
#      key = int(re.sub(r'f([0-9].*)', '\\1', keys[i]))
#      f_i = key // feature_vecsize
#      f_ij = key % feature_vecsize
#      keys_names.append('%s_%d' % (features[f_i], f_ij))
#
#  plt.figure(figsize=(10, 6))
#  plt.title('Feature Importance by Gain')
#  plt.bar(keys_names, [values[i] for i in indices])
#  plt.xticks(rotation=90)
#  plt.tight_layout()
#  plt.savefig(outfolder+'feature_importance.pdf')
#
#  print("Importances : ")
#  for i in indices:
#      print(keys_names[i], values[i])
