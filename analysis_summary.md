# Analysis Scripts Summary

## `examples/c++/pmt_bdt.cc` — Signal / Physics Events

Reads RAT::DS Monte Carlo output (muon events) and produces a DOM-level feature tree for BDT training.

**Inputs:**
- RAT ROOT files matching a glob pattern
- `xsec_graphs.root` — neutrino cross-section graphs for rock-bed event re-weighting

**Processing:**
1. Reads PMT geometry (positions) from the `runT` tree to compute mean DOM (x, y, z).
2. Loops over triggered events; extracts muon vertex, momentum, and kinetic energy from MC truth.
3. Skips type-0 (8-inch HQE) PMTs; keeps type-1 (12-inch) only.
4. Groups hit PMTs into DOMs (31 PMTs per DOM: `dom_id = pmt_id / 31`).
5. Per DOM: sorts PE and time vectors, runs a rolling consecutive-hit window to compute `npmts_50/100` and `npe_50/100`, then computes PE and timing statistics.

**Output TTree `doms`** (one entry per DOM per event) + `TH1D nevts` (weighted event count):

| Branch group | Branches |
|---|---|
| Event truth | `event_id`, `vtxX/Y/Z`, `momX/Y/Z`, `muE`, `rock_wgt` |
| Geometry | `dom_id`, `du_id`, `dom_x/y/z` |
| PE stats | `npmts`, `npe`, `pe_min`, `pe_spread`, `pe_rms` |
| Timing stats | `t_min`, `t_spread`, `t_mean`, `t_rms` |
| Coincidence | `npmts_50`, `npe_50`, `npmts_100`, `npe_100` |

**Build & run:**
```bash
cd examples/c++
make pmt_bdt
./pmt_bdt "/path/to/ratfiles/*.root" output_signal.root
```

---

## `examples/c++/ambient_bdt.cc` — Ambient Light Background

Generates synthetic ambient-light background events by randomly sampling per-PMT time frames from a pre-simulated `frame_light` TTree. No RAT dependency — ROOT only.

**Input:**
- A ROOT file containing a `frame_light` TTree (e.g. `output_ambient_light_10us.root`)
  - Each entry holds one PMT's hits in a fixed-length time window
  - Only `mcpehittime_rel` (hit time relative to frame start, ns) is read

**Processing:**
1. **Pre-loads** all 5M frame entries into a compact `vector<FrameData>` (`npe`, `t_first` per entry). ~40 MB in memory; no PMT ID comparison needed.
2. For each of **M iterations** (synthetic events):
   - Draws a random frame index `seq[i] ∈ [0, nentries)` independently for each of the 2790 PMTs.
   - Groups sampled hits into DOMs and runs the **identical rolling-window and statistics analysis** as `pmt_bdt.cc`.
3. Fills the output `doms` tree with the same branch structure as `pmt_bdt.cc` (muon-truth branches omitted; `dom_x/y/z = 0`).

**Output TTree `doms`** (one entry per hit DOM per event):

| Branch group | Branches |
|---|---|
| Event index | `event_id` |
| Geometry | `dom_id`, `du_id`, `dom_x/y/z` (zero) |
| PE stats | `npmts`, `npe`, `pe_min`, `pe_spread`, `pe_rms` |
| Timing stats | `t_min`, `t_mean`, `t_rms` |
| Coincidence | `npmts_50`, `npe_50`, `npmts_100`, `npe_100` |

**Build & run:**
```bash
cd examples/c++
make ambient_bdt
./ambient_bdt <input_root> <output_root> <M> [seed]

# Example — 1000 events, time-based random seed:
./ambient_bdt ../../output_ambient_light_10us.root ambient_out.root 1000

# Example — reproducible run with fixed seed:
./ambient_bdt ../../output_ambient_light_10us.root ambient_out.root 1000 42
```

---

## Output Compatibility

Both scripts write a TTree named `doms` with the same observable branches, so signal and background outputs can be combined directly for BDT training:

```bash
# In ROOT:
TChain ch("doms");
ch.Add("output_signal.root");
ch.Add("ambient_out.root");
```
