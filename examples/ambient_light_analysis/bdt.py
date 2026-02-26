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

frates = ROOT.TFile("/home/nitish/public_html/shared/uboone/flugg_study/for_milind/plots_pretty/data/event_rates_lp3_depth_numu.root", "read")
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

fsignal = ROOT.TFile("bdtinput_signal_aframe_spacing5m_hex_DU_v2.root", "read")
fbkg = ROOT.TFile("bdtinput_bkg_aframe_spacing5m_hex_DU_v2.root", "read")
tot_signal = fsignal.Get("nevts").Integral()
tot_bkg = fbkg.Get("nevts").Integral()

# now try to normalize everything to per Day
pot = tot_signal/rate_perPOT
osc_factor = 0.73
pot_perDay = 7.92E-3
norm_signal = osc_factor*pot_perDay/pot

livetime = tot_bkg/cosmic_rate
spill_rate = 0.937 # Hz
spill_length = 1.E-5 # 10us
norm_bkg = 86400*spill_rate*spill_length/livetime

print("Normalization Factors : ", norm_signal, norm_bkg)

outfolder = "/home/nitish/public_html/shared/end/efficiency/domPE/bdt_aframe_spacing5m_hex_DU_v2/3dom_3pe_50/"

# Load dataset
X_signal = np.load('signal_aframe_spacing5m_hex_DU_v2.npy')
X_bkg = np.load('bkg_aframe_spacing5m_hex_DU_v2.npy')
X_signal[np.isinf(X_signal)] = -5
X_bkg[np.isinf(X_bkg)] = -5

# aux info
df_signal = pd.read_csv('signal_aframe_spacing5m_hex_DU_v2_info.csv').to_numpy()
df_bkg = pd.read_csv('bkg_aframe_spacing5m_hex_DU_v2_info.csv').to_numpy()
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
X = np.vstack([X_signal, X_bkg])              # shape (n_signal + n_bkg, L * n_cols)
y = np.concatenate([y_signal, y_bkg])         # shape (n_signal + n_bkg,)
del X_signal, X_bkg

# Split: 60% train, 20% validation, 20% test
X_train, X_temp, y_train, y_temp, info_train, info_temp = train_test_split(
    X, y, X_info, test_size=0.2, stratify=y, random_state=42
)
X_val, X_test, y_val, y_test, info_val, info_test = train_test_split(
    X_temp, y_temp, info_temp, test_size=0.5, stratify=y_temp, random_state=43
)

#  signal_evtrate_factor = X_info[:,0][y == 1].sum()/info_test[:,0][y_test == 1].sum()
#  bkg_evtrate_factor = X_info[:,0][y == 0].sum()/info_test[:,0][y_test == 0].sum()
#  print("Test vs Total Factors : ", signal_evtrate_factor, bkg_evtrate_factor)
signal_evtrate_factor = 1.
bkg_evtrate_factor = 1.


# Create DMatrix objects
dtrain = xgb.DMatrix(X_train, label=y_train)
dval = xgb.DMatrix(X_val, label=y_val)
#  dtest = xgb.DMatrix(X_test, label=y_test)
# For stats predict on full set
# we've checked performance is not too different from training and test
dtest = xgb.DMatrix(X, label=y)
y_test = y
info_test = X_info

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
    'scale_pos_weight': 3
}

# Train with early stopping on validation set
evals = [(dtrain, 'train'), (dval, 'val')]
start=time.time()
bst = xgb.train(
    params,
    dtrain,
    num_boost_round=500,
    evals=evals,
    early_stopping_rounds=20,
    verbose_eval=10
)
print("Time Taken for Training (s) : ", time.time()-start)

# Save model for later inference
bst.save_model('xgb_model_aframe_spacing5m_hex_DU_v2.json')

def apply_trigger(X):
    #  # older tree
    #  ndoms_pecut = np.sum(X[:, 3*feature_vecsize:4*feature_vecsize] >= 3, axis=1)
    #  # newer tree, 100ns for pe
    #  ndoms_pecut = np.sum(X[:, 6*feature_vecsize:7*feature_vecsize] >= 3, axis=1)
    # newer tree, 50ns for pe
    ndoms_pecut = np.sum(X[:, 7*feature_vecsize:8*feature_vecsize] >= 3, axis=1)
    return np.where(ndoms_pecut >=3)

# Evaluate on the test set
y_pred_prob = bst.predict(dtest)
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

trigger_cut = apply_trigger(X)
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
