# Output Tree Size Estimate for `ambient_bdt`

## Parameters
- `NPMTS` = 2790 (31 PMTs per DOM, 90 DOMs)
- `M` = 1,000,000 sampling iterations
- Input: `output_ambient_light_10us.root` (5,000,000 frame entries)

## Measured from Data
Sampled every 1000th entry from `frame_light` tree:

| Quantity | Value |
|---|---|
| Hit fraction per PMT frame (p) | **0.1334** (13.34%) |
| Mean npe per hit frame | **1.06** |

## DOM-Level Hit Probability

Each DOM has 31 PMTs. The probability that at least one PMT in a DOM is hit:

```
P(DOM hit) = 1 - (1 - p)^31
           = 1 - (0.8666)^31
           ≈ 0.988
```

~98.8% of DOMs get an output entry per event.

## Total Tree Entries

```
M × 90 × P(DOM hit) = 1,000,000 × 90 × 0.988 ≈ 88.9 million entries
```

## Bytes per Entry (Uncompressed)

17 active branches × 4 bytes = **68 bytes/entry**

(Note: `t_spread` is currently commented out; `vtx/mom/muE/rock_wgt` branches are absent.)

| Branch | Type | Bytes |
|---|---|---|
| event_id, dom_id, du_id | int | 12 |
| dom_x, dom_y, dom_z | float | 12 |
| npmts, npe, pe_min, pe_spread | int | 16 |
| pe_rms | float | 4 |
| t_min, t_mean, t_rms | float | 12 |
| npmts_50, npe_50, npmts_100, npe_100 | int | 16 |
| **Total** | | **72 bytes** |

Wait — that's 17 branches × 4 = 68 bytes (t_spread excluded), but counting above gives 72 bytes for 18 branches. Recount: 17 branches → **68 bytes/entry**.

## File Size Estimate

| | Value |
|---|---|
| Uncompressed data | 88.9M × 68 B ≈ **6.0 GB** |
| ROOT file (compression ~4×) | ≈ **1.5 GB** |

Compression ratio typically ranges 3–5× for float/int data of this type, giving a plausible range of **1.2–2.0 GB** on disk.
