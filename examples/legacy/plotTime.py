import uproot
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

def plot_pmt_event(filename, event_index=0, projection='XY',
                   scale_factor=50, base_size=10, cmap='viridis',
                   tmin=None, tmax=None, no_charge_scaling=False):
    """
    Plots PMT hits for a specific event from a ratpac-two ntuple file,
    with optional time filtering and charge scaling.

    Args:
        filename (str): Path to the ROOT ntuple file.
        event_index (int): Index of the event to plot (0-based).
        projection (str): The 2D plane to project onto ('XY', 'YZ', 'XZ').
        scale_factor (float): Multiplier for charge to determine marker size increase
                              (used only if no_charge_scaling is False).
        base_size (float): Base marker size (used as minimum if scaling, or uniform size if not).
        cmap (str): Matplotlib colormap name for time visualization.
        tmin (float, optional): Minimum hit time to plot. Defaults to None (no minimum).
        tmax (float, optional): Maximum hit time to plot. Defaults to None (no maximum).
        no_charge_scaling (bool): If True, disable marker size scaling by charge.
                                  Defaults to False.
    """
    try:
        file = uproot.open(filename)
    except Exception as e:
        print(f"Error opening file {filename}: {e}")
        sys.exit(1)

    # --- 1. Read Meta Tree for PMT Positions ---
    if "meta" not in file:
        print(f"Error: 'meta' TTree not found in {filename}")
        sys.exit(1)
    meta_tree = file["meta"]

    if meta_tree.num_entries == 0:
        print(f"Error: 'meta' TTree in {filename} has no entries.")
        sys.exit(1)

    try:
        meta_data = meta_tree.arrays(
            ["pmtId", "pmtX", "pmtY", "pmtZ"],
            library="np",
            entry_stop=1
        )
        pmt_ids_meta = meta_data["pmtId"][0]
        pmt_x = meta_data["pmtX"][0]
        pmt_y = meta_data["pmtY"][0]
        pmt_z = meta_data["pmtZ"][0]
        pmt_locations = {pid: (x, y, z) for pid, x, y, z in zip(pmt_ids_meta, pmt_x, pmt_y, pmt_z)}
        print(f"Read location data for {len(pmt_locations)} PMTs from 'meta' tree.")
    except KeyError as e:
        print(f"Error: Missing branch in 'meta' TTree: {e}")
        print("Required branches: pmtId, pmtX, pmtY, pmtZ")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading 'meta' tree: {e}")
        sys.exit(1)

    # --- 2. Read Output Tree for Event Data ---
    if "output" not in file:
        print(f"Error: 'output' TTree not found in {filename}")
        sys.exit(1)
    output_tree = file["output"]

    if event_index >= output_tree.num_entries:
        print(f"Error: Event index {event_index} is out of range ({output_tree.num_entries} events).")
        sys.exit(1)

    try:
        event_data = output_tree.arrays(
            ['hitPMTID', 'hitPMTTime', 'hitPMTCharge'],
            entry_start=event_index,
            entry_stop=event_index + 1,
            library="np"
        )
        if len(event_data['hitPMTID']) == 0:
             print(f"Error: Could not read data for event index {event_index}.")
             sys.exit(1)

        hit_ids = event_data['hitPMTID'][0]
        hit_times = event_data['hitPMTTime'][0]
        hit_charges = event_data['hitPMTCharge'][0]
        print(f"Read {len(hit_ids)} total hits for event {event_index} from 'output' tree.")

    except KeyError as e:
        print(f"Error: Missing branch in 'output' TTree: {e}")
        print("Required branches: hitPMTID, hitPMTTime, hitPMTCharge")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading event {event_index} from 'output' tree: {e}")
        sys.exit(1)

    # --- 3. Apply Time Filtering (if requested) ---
    time_mask = np.full(len(hit_times), True) # Start with all True
    filter_applied = False
    if tmin is not None:
        time_mask &= (hit_times >= tmin)
        print(f"Applying time filter: t >= {tmin}")
        filter_applied = True
    if tmax is not None:
        time_mask &= (hit_times <= tmax)
        print(f"Applying time filter: t <= {tmax}")
        filter_applied = True

    if filter_applied:
        original_hit_count = len(hit_ids)
        hit_ids = hit_ids[time_mask]
        hit_times = hit_times[time_mask]
        hit_charges = hit_charges[time_mask]
        print(f"Filtered hits by time: {original_hit_count} -> {len(hit_ids)}")

    if len(hit_ids) == 0:
        print(f"Event {event_index} has no recorded hits within the specified criteria. Nothing to plot.")
        sys.exit(0)

    # --- 4. Prepare Data for Plotting ---
    plot_coords1 = []
    plot_coords2 = []
    plot_times = []
    plot_charges = []
    missing_pmts = 0

    projection = projection.upper()
    if projection == 'XY':
        idx1, idx2 = 0, 1
        axis_label1, axis_label2 = 'X', 'Y'
    elif projection == 'YZ':
        idx1, idx2 = 1, 2
        axis_label1, axis_label2 = 'Y', 'Z'
    elif projection == 'XZ':
        idx1, idx2 = 0, 2
        axis_label1, axis_label2 = 'X', 'Z'
    else:
        print(f"Error: Invalid projection '{projection}'. Choose from 'XY', 'YZ', 'XZ'.")
        sys.exit(1)

    # Map hit IDs to locations and collect data for valid, filtered hits
    for i, pmt_id in enumerate(hit_ids):
        if pmt_id in pmt_locations:
            coords = pmt_locations[pmt_id]
            plot_coords1.append(coords[idx1])
            plot_coords2.append(coords[idx2])
            plot_times.append(hit_times[i])
            plot_charges.append(hit_charges[i])
        else:
            missing_pmts += 1

    if missing_pmts > 0:
        print(f"Warning: {missing_pmts} hit PMT IDs (after filtering) were not found in the meta data and are excluded.")

    if not plot_coords1:
        print("No valid PMT hits found to plot after matching with meta data and filtering.")
        sys.exit(0)

    plot_coords1 = np.array(plot_coords1)
    plot_coords2 = np.array(plot_coords2)
    plot_times = np.array(plot_times)
    plot_charges = np.array(plot_charges)

    # --- 5. Calculate Marker Sizes ---
    plot_charges = np.maximum(plot_charges, 0) # Ensure non-negative charges

    if no_charge_scaling:
        print("Charge scaling disabled. Using uniform marker size.")
        marker_sizes = np.full_like(plot_charges, base_size)
    else:
        print(f"Scaling marker area by charge (base: {base_size}, scale: {scale_factor}).")
        # s parameter is proportional to area. Linear scaling with charge:
        marker_sizes = base_size + plot_charges * scale_factor

    # --- 6. Create Plot ---
    fig, ax = plt.subplots(figsize=(10, 8))

    # Determine vmin and vmax for color scale normalization
    # If user specifies tmin/tmax, use those for the color scale limits.
    # Otherwise, let matplotlib determine from the plotted data range.
    vmin_plot = tmin if tmin is not None else np.min(plot_times) if len(plot_times) > 0 else 0
    vmax_plot = tmax if tmax is not None else np.max(plot_times) if len(plot_times) > 0 else 1
    # Handle edge case where filter results in single time value or tmin == tmax
    if vmin_plot == vmax_plot:
        print("Warning: Time range for color scale has zero width. Adjusting slightly.")
        vmin_plot -= 0.5
        vmax_plot += 0.5

    scatter = ax.scatter(
        plot_coords1,
        plot_coords2,
        s=marker_sizes,
        c=plot_times,
        cmap=cmap,
        alpha=0.7,
        edgecolors='k',
        linewidths=0.5,
        vmin=vmin_plot, # Set color scale minimum
        vmax=vmax_plot  # Set color scale maximum
    )

    cbar = fig.colorbar(scatter)
    time_unit = "ns" # ASSUMPTION: Modify if your units are different
    cbar.set_label(f'Hit Time ({time_unit})')

    pos_unit = "mm" # ASSUMPTION: Modify if your units are different
    ax.set_xlabel(f'{axis_label1} Coordinate ({pos_unit})')
    ax.set_ylabel(f'{axis_label2} Coordinate ({pos_unit})')
    title_str = f'PMT Hits - Event {event_index} ({projection} Projection)'
    if tmin is not None or tmax is not None:
        t_range_str = f"Time Range: [{tmin if tmin is not None else '-inf'}, {tmax if tmax is not None else 'inf'}] {time_unit}"
        title_str += f'\n{t_range_str}'
    ax.set_title(title_str) # Removed filename from title for brevity

    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, linestyle='--', alpha=0.6)

    num_hits_plotted = len(plot_coords1)
    ax.text(0.02, 0.02, f'Plotted Hits: {num_hits_plotted}', transform=ax.transAxes,
            fontsize=9, verticalalignment='bottom', bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.5))

    print(f"Displaying plot for event {event_index}, projection {projection}.")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot PMT hits from a ratpac-two ROOT ntuple file.")
    parser.add_argument("filename", help="Path to the input ROOT ntuple file.")
    parser.add_argument("-e", "--event", type=int, default=0, dest="event_index",
                        help="Index of the event to plot (default: 0).")
    parser.add_argument("-p", "--projection", type=str, default="XY", choices=["XY", "YZ", "XZ"],
                        help="2D projection plane (default: XY).")
    parser.add_argument("--scale", type=float, default=50, dest="scale_factor",
                        help="Scale factor for charge -> marker area increase (default: 50). Ignored if --no_charge_scaling is set.")
    parser.add_argument("--base_size", type=float, default=10,
                        help="Base marker size (default: 10). Used as minimum or uniform size.")
    parser.add_argument("--cmap", type=str, default="viridis",
                        help="Matplotlib colormap for time (default: viridis). E.g., 'plasma', 'jet', 'coolwarm'.")
    parser.add_argument("--tmin", type=float, default=None,
                        help="Minimum hit time to plot (default: no minimum).")
    parser.add_argument("--tmax", type=float, default=None,
                        help="Maximum hit time to plot (default: no maximum).")
    parser.add_argument("--no_charge_scaling", action="store_true",
                        help="If set, disable marker area scaling by charge (use uniform base_size).")

    args = parser.parse_args()

    # Optional: Add validation if tmin and tmax are both specified
    if args.tmin is not None and args.tmax is not None and args.tmin > args.tmax:
        print(f"Warning: tmin ({args.tmin}) is greater than tmax ({args.tmax}). No hits will be selected.")
        # Or could exit: sys.exit("Error: tmin cannot be greater than tmax.")

    plot_pmt_event(
        filename=args.filename,
        event_index=args.event_index,
        projection=args.projection,
        scale_factor=args.scale_factor,
        base_size=args.base_size,
        cmap=args.cmap,
        tmin=args.tmin,
        tmax=args.tmax,
        no_charge_scaling=args.no_charge_scaling
    )
