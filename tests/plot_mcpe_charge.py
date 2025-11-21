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
# # Visualize MC photo-electrons
#
# Extracts Monte Carlo PE positions and charges from the `output` tree in
# `muon_balabala.root` and renders them as a 3D scatter plot with colors mapped
# to `mcPECharge`.

# %%
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

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
        description="Visualize Monte Carlo PE positions colored by charge."
    )
    parser.add_argument(
        "--input-file",
        type=Path,
        default=Path("muon_balabala.root"),
        help="ROOT file containing the 'output' tree.",
    )
    parser.add_argument(
        "--event-index",
        type=int,
        default=0,
        help="Zero-based event index to visualize.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=0.3,
        help="Marker size for the k3d scatter points.",
    )
    parser.add_argument(
        "--min-charge",
        type=float,
        default=None,
        help="Optional lower threshold on mcPECharge to keep hits.",
    )
    # Ignore spurious notebook arguments (e.g., -f <connection.json>).
    args, _ = parser.parse_known_args(cli_args)
    return args


def fetch_event(
    output_tree: uproot.behaviors.TBranch.HasBranches, event_index: int
) -> dict[str, ak.Array]:
    total_entries = output_tree.num_entries
    if event_index < 0 or event_index >= total_entries:
        raise IndexError(
            f"Requested event {event_index} but tree holds {total_entries} entries"
        )
    return output_tree.arrays(
        ["mcPEPMTID", "mcPEx", "mcPEy", "mcPEz", "mcPECharge"],
        entry_start=event_index,
        entry_stop=event_index + 1,
        library="ak",
    )


def flatten_event(event_arrays) -> dict[str, np.ndarray]:
    """Convert Awkward record-of-arrays into flat numpy arrays per field."""

    return {
        field: ak.to_numpy(ak.flatten(event_arrays[field]))
        for field in event_arrays.fields
    }


def apply_charge_cut(
    data: dict[str, np.ndarray], min_charge: float | None
) -> dict[str, np.ndarray]:
    if min_charge is None:
        return data
    mask = data["mcPECharge"] >= float(min_charge)
    filtered = {key: values[mask] for key, values in data.items()}
    return filtered


def make_colors(values: np.ndarray) -> tuple[np.ndarray, tuple[float, float]]:
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
        rgba = get_cmap("inferno")(norm)
        rgb = np.round(rgba[:, :3] * 255).astype(np.uint8)
    else:
        rgb = np.stack([norm, 1.0 - norm, 0.5 + 0.5 * norm], axis=1)
        rgb = np.round(np.clip(rgb, 0.0, 1.0) * 255).astype(np.uint8)
    colors = (
        (rgb[:, 0].astype(np.uint32) << 16)
        | (rgb[:, 1].astype(np.uint32) << 8)
        | rgb[:, 2].astype(np.uint32)
    )
    return colors, (vmin, vmax)


def visualize_mcpe(
    file_path: Path,
    event_index: int = 0,
    point_size: float = 0.3,
    min_charge: float | None = None,
) -> k3d.Plot:
    file_path = Path(file_path).expanduser().resolve()
    if not file_path.exists():
        raise FileNotFoundError(file_path)

    with uproot.open(file_path) as infile:
        try:
            output_tree = infile["output"]
        except KeyError as exc:
            raise KeyError("Input file lacks 'output' tree") from exc
        event_arrays = fetch_event(output_tree, event_index)

    flat_event = flatten_event(event_arrays)
    filtered = apply_charge_cut(flat_event, min_charge)
    coords = np.column_stack(
        [filtered["mcPEx"], filtered["mcPEy"], filtered["mcPEz"]]
    ).astype(np.float32)
    charges = filtered["mcPECharge"].astype(np.float32)

    if coords.size == 0:
        raise RuntimeError("No MC PE coordinates left after filtering.")

    colors, color_range = make_colors(charges)
    plot = k3d.plot(name=f"MC PE charges event {event_index}")
    plot += k3d.points(
        coords,
        point_size=point_size,
        colors=colors,
        shader="3d",
    )
    plot.axes = ["X (cm)", "Y (cm)", "Z (cm)"]

    vmin, vmax = color_range
    if np.isfinite(vmin) and np.isfinite(vmax):
        print(f"mcPECharge range: {vmin:.3f} -- {vmax:.3f}")
    if min_charge is not None:
        print(
            f"Applied charge cut at {min_charge}; kept {charges.size} hits from"
            f" {flat_event['mcPECharge'].size}."
        )

    return plot


def main(cli_args: Iterable[str] | None = None) -> k3d.Plot:
    args = parse_args(cli_args)
    plot = visualize_mcpe(
        file_path=args.input_file,
        event_index=args.event_index,
        point_size=args.point_size,
        min_charge=args.min_charge,
    )
    plot.display()
    return plot


if __name__ == "__main__":
    main()
