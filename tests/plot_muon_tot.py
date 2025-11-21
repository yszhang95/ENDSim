#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Visualize `digitTimeOverThreshold` for the first muon event
#
# Loads PMT geometry from the `meta` tree and overlays the first event in
# `muon_balabala.root` as a colored 3D point cloud using `k3d`.

# %%
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, Tuple

import awkward as ak
import numpy as np
import k3d
import uproot

try:
    from matplotlib import colormaps as mpl_colormaps

    def get_cmap(name: str):
        return mpl_colormaps[name]

except Exception:
    try:
        from matplotlib import cm as mpl_cm

        def get_cmap(name: str):
            return mpl_cm.get_cmap(name)

    except Exception:  # matplotlib is optional; fall back to a simple gradient
        get_cmap = None


def parse_args(cli_args: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a 3D view of digitized PMT hits for a single event using k3d."
    )
    parser.add_argument(
        "--input-file",
        type=Path,
        default=Path("muon_1gev_output.root"),
        help="Path to the ROOT file containing 'output' and 'meta' trees.",
    )
    parser.add_argument(
        "--event-index",
        type=int,
        default=0,
        help="Zero-based event index to visualize (default: 0, the first event).",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=0.4,
        help="Marker size passed to k3d.points (default: 0.4).",
    )
    # Jupyter injects its own -f/--f arguments when running a cell; ignore them.
    args, _ = parser.parse_known_args(cli_args)
    return args


def load_geometry(meta_tree: uproot.behaviors.TBranch.HasBranches) -> Dict[int, np.ndarray]:
    """Return a mapping of PMT id -> (x, y, z)."""

    meta_arrays = meta_tree.arrays(["pmtId", "pmtX", "pmtY", "pmtZ"], library="ak")
    pmt_ids = ak.to_numpy(ak.flatten(meta_arrays["pmtId"]))
    pmt_x = ak.to_numpy(ak.flatten(meta_arrays["pmtX"]))
    pmt_y = ak.to_numpy(ak.flatten(meta_arrays["pmtY"]))
    pmt_z = ak.to_numpy(ak.flatten(meta_arrays["pmtZ"]))
    coords = np.column_stack([pmt_x, pmt_y, pmt_z]).astype(np.float32)
    return {int(pid): coords[idx] for idx, pid in enumerate(pmt_ids)}


def fetch_event(
    output_tree: uproot.behaviors.TBranch.HasBranches, event_index: int
) -> dict[str, ak.Array]:
    total_entries = output_tree.num_entries
    if event_index < 0 or event_index >= total_entries:
        raise IndexError(
            f"Requested event {event_index} but tree holds {total_entries} entries"
        )
    return output_tree.arrays(
        ["digitPMTID", "digitTimeOverThreshold"],
        entry_start=event_index,
        entry_stop=event_index + 1,
        library="ak",
    )


def attach_positions(
    digit_ids: np.ndarray, tot_values: np.ndarray, geometry: Dict[int, np.ndarray]
) -> Tuple[np.ndarray, np.ndarray, set[int]]:
    coords: list[np.ndarray] = []
    kept_tot: list[float] = []
    missing: set[int] = set()
    for pid, tot in zip(digit_ids, tot_values):
        pos = geometry.get(int(pid))
        if pos is None:
            missing.add(int(pid))
            continue
        coords.append(pos)
        kept_tot.append(float(tot))
    if coords:
        coord_array = np.vstack(coords).astype(np.float32)
        tot_array = np.asarray(kept_tot, dtype=np.float32)
    else:
        coord_array = np.empty((0, 3), dtype=np.float32)
        tot_array = np.empty((0,), dtype=np.float32)
    return coord_array, tot_array, missing


def make_colors(values: np.ndarray) -> Tuple[np.ndarray, Tuple[float, float]]:
    if values.size == 0:
        return np.empty((0,), dtype=np.uint32), (np.nan, np.nan)
    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    span = vmax - vmin
    if not np.isfinite(span) or span <= 0.0:
        norm = np.zeros_like(values, dtype=np.float32)
    else:
        norm = (values - vmin) / span
    if get_cmap is not None:
        rgba = get_cmap("viridis")(norm)
        rgb = np.round(rgba[:, :3] * 255).astype(np.uint8)
    else:
        rgb = np.stack([norm, 0.5 * (1.0 - np.abs(norm - 0.5) * 2), 1.0 - norm], axis=1)
        rgb = np.round(np.clip(rgb, 0.0, 1.0) * 255).astype(np.uint8)
    colors = (
        (rgb[:, 0].astype(np.uint32) << 16)
        | (rgb[:, 1].astype(np.uint32) << 8)
        | rgb[:, 2].astype(np.uint32)
    )
    return colors, (vmin, vmax)


def visualize_event(
    file_path: Path, event_index: int = 0, point_size: float = 0.4
) -> k3d.Plot:
    file_path = Path(file_path).expanduser().resolve()
    if not file_path.exists():
        raise FileNotFoundError(file_path)

    with uproot.open(file_path) as infile:
        try:
            output_tree = infile["output"]
            meta_tree = infile["meta"]
        except KeyError as exc:
            raise KeyError("Input file must contain 'output' and 'meta' trees") from exc

        geometry = load_geometry(meta_tree)
        event = fetch_event(output_tree, event_index)

    digit_ids = ak.to_numpy(ak.flatten(event["digitPMTID"]))
    tot_values = ak.to_numpy(ak.flatten(event["digitTimeOverThreshold"]))
    coords, filtered_tot, missing = attach_positions(digit_ids, tot_values, geometry)
    print(f'Visualizing event {event_index} from file: {file_path}')
    print(f'Total PMT hits: {len(digit_ids)}, Plotted hits: {coords.shape[0]}')
    if coords.size == 0:
        raise RuntimeError(
            "No PMT coordinates found for the requested event; check geometry contents."
        )

    colors, color_range = make_colors(filtered_tot)
    plot = k3d.plot(name=f"digitTimeOverThreshold event {event_index}")
    plot += k3d.points(
        coords,
        point_size=point_size,
        colors=colors,
        shader="3d",
    )
    plot.axes = ["X (cm)", "Y (cm)", "Z (cm)"]

    vmin, vmax = color_range
    if np.isfinite(vmin) and np.isfinite(vmax):
        print(f"digitTimeOverThreshold range: {vmin:.3f} -- {vmax:.3f}")
    if missing:
        print(f"Skipped {len(missing)} PMT ids without geometry: {sorted(missing)[:5]}...")

    return plot


def main(cli_args: Iterable[str] | None = None) -> k3d.Plot:
    args = parse_args(cli_args)
    plot = visualize_event(args.input_file, args.event_index, args.point_size)
    plot.display()
    return plot


if __name__ == "__main__":
    main()
