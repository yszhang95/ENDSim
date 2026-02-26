# Changes to `bdt.py` — Streaming Data Loading & Dynamic scale_pos_weight

## Overview

`bdt.py` was modified to:
1. Read event data directly from ROOT files at training time (no intermediate `.npy`/`.csv` on disk).
2. Limit how many signal and background events are loaded via CLI arguments.
3. Compute `scale_pos_weight` automatically from the loaded sample sizes, or accept an explicit override.

All other training logic (normalization, trigger, train/val/test split, XGBoost params, plots) is unchanged.

---

## Pipeline Context

```
RAT-PAC simulation (.root files)
         |
         v
  pmt_bdt.cc   (C++)
  --> doms TTree: one row per active DOM per event
  --> nevts histogram: total weighted event count
         |
         v
  bdt.py  (previously needed save_csv.py first)
  --> now reads ROOT directly via stream_events_from_root()
  --> xgb_model_*.json, ROC/efficiency/purity plots
```

`save_csv.py` is no longer a required preprocessing step. It can still be used independently if persistent `.npy` files are needed.

---

## New CLI Arguments

```
python bdt.py \
  --signal_root  <path>   # required: signal ROOT file from pmt_bdt.cc
  --bkg_root     <path>   # required: background ROOT file from pmt_bdt.cc
  --n_signal     <int>    # required: max signal events to load
  --n_bkg        <int>    # required: max background events to load
  --scale_pos_weight <float>  # optional: XGBoost weight (default: n_bkg/n_signal)
  --max_len      <int>    # optional: DOM padding slots per event (default: 88)
```

### Example

```bash
cd examples/ambient_light_analysis
python bdt.py \
  --signal_root ../c++/bdtinput_signal_numu.root \
  --bkg_root    ../c++/bdtinput_bkg.root \
  --n_signal 5000 \
  --n_bkg 15000
# scale_pos_weight auto-computed as 15000/5000 = 3.0

# Override:
python bdt.py ... --scale_pos_weight 5.0
```

---

## New Function: `stream_events_from_root`

```python
def stream_events_from_root(root_file, max_events, max_len=88):
```

**What it does:**

1. Opens `root_file` and reads the entire `doms` TTree into memory using `ROOT.RDataFrame(...).AsNumpy()` — vectorized C++ reads, much faster than a Python `for entry in tree:` loop.
2. Groups all rows by `event_id` using `pandas.groupby`:
   - Info columns (`rock_wgt`, `vtxX/Y/Z`, `momX/Y/Z`, `muE`): averaged across DOM rows (they are identical for all DOMs in an event).
   - Feature columns (all other branches): collected as a list of per-DOM values.
3. Takes the first `max_events` grouped events via `iloc[:max_events]`.
4. Pads or truncates each feature column list to `max_len` slots (padding value: `-1.0`).
5. Stacks all feature columns horizontally into a `(n_events, n_feat_cols * max_len)` float32 array.

**Returns:**
- `X`: `np.float32` array, shape `(n_events, n_feat_cols * max_len)`
- `info`: `np.float32` array, shape `(n_events, 8)`, columns `[rock_wgt, vtxX, vtxY, vtxZ, momX, momY, momZ, muE]`

**Feature layout** is determined by the TTree branch order (preserved by RDataFrame), so the feature index arithmetic in `apply_trigger` (`X[:, 7*feature_vecsize : 8*feature_vecsize]`) remains valid as long as `max_len == feature_vecsize == 88`.

---

## `scale_pos_weight` Logic

```python
if args.scale_pos_weight is not None:
    scale_pos_weight = args.scale_pos_weight
else:
    scale_pos_weight = len(y_bkg) / max(len(y_signal), 1)
```

Previously hardcoded as `3`. Now derived from the actual loaded sample sizes so it stays correct if `--n_signal` / `--n_bkg` are changed.

---

## What Was Removed

- Loading from `signal_*.npy` / `bkg_*.npy`
- Loading from `signal_*_info.csv` / `bkg_*_info.csv`
- Hardcoded `bdtinput_signal_*.root` / `bdtinput_bkg_*.root` paths (replaced by `args.signal_root` / `args.bkg_root`)
- Hardcoded `scale_pos_weight: 3` in XGBoost params

---

## Unchanged

- All normalization constants (`norm_signal`, `norm_bkg`, `rate_perPOT`, `cosmic_rate`, etc.)
- Hardcoded event-rate ROOT file path (`/home/nitish/.../event_rates_lp3_depth_numu.root`)
- `apply_trigger()`, `angle()` functions
- 60/20/20 stratified train/val/test split
- XGBoost training loop, early stopping, model save
- All diagnostic plots (ROC, BDT score vs angle/energy, FOM optimization)
- `outfolder` for plot output
