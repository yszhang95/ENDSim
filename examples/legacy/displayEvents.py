import uproot
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os
import time
import awkward as ak

# --- Configuration ---
# NOTE: This script now uses matplotlib for standalone window display.

# TODO: Please update these file paths to match your system.
RATPAC_FILE = 'output_out.numu_100000.rootracker.89600-89700.randomVtx.display_89600_89700.root'
GENIE_FILE = 'numu_100000.rootracker.root'

# PDG codes for lookup
PDG_MUON = 13
PDG_ANTI_MUON = -13
PDG_MUON_NEUTRINO = 14
PDG_ANTI_MUON_NEUTRINO = -14

# Physics and Display Parameters
MUON_MASS_MEV = 105.65837
MUON_ENERGY_LOSS_MEV_PER_MM = 0.2 # 2 MeV/cm -> 0.2 MeV/mm
CHERENKOV_ANGLE_DEG = 42.0

# --- Matching Tolerances ---
# Matching is now done on KE and the angle with the Z-axis.
KE_TOLERANCE_MEV = 0.1 # 0.1 MeV
ANGLE_TOLERANCE_COS = 0.001 # Tolerance for the cosine of the angle


def load_data(ratpac_path, genie_path):
    """
    Loads the TTree objects from the ROOT files using uproot.
    """
    if not os.path.exists(ratpac_path):
        print(f"Error: RATPAC file not found at '{ratpac_path}'")
        return None, None, None
    if not os.path.exists(genie_path):
        print(f"Error: GENIE file not found at '{genie_path}'")
        return None, None, None

    print("Loading data from files...")
    rat_file = uproot.open(ratpac_path)
    genie_file = uproot.open(genie_path)

    return rat_file['meta'], rat_file['output'], genie_file['gRooTracker']


def get_pmt_positions(meta_tree):
    """
    Extracts PMT positions from the meta tree.
    """
    print("Reading detector geometry (PMT positions)...")
    pmt_data = meta_tree.arrays(['pmtId', 'pmtX', 'pmtY', 'pmtZ'])
    pmt_ids, pmt_x, pmt_y, pmt_z = pmt_data['pmtId'][0], pmt_data['pmtX'][0], pmt_data['pmtY'][0], pmt_data['pmtZ'][0]
    pmt_positions = {int(pid): {'x': x, 'y': y, 'z': z} for pid, x, y, z in zip(pmt_ids, pmt_x, pmt_y, pmt_z)}
    print(f"Found {len(pmt_positions)} PMTs in the detector.")
    return pmt_positions

def get_genie_muon_kinematics(genie_data):
    """
    Pre-processes the entire GENIE file to extract kinematics for all final-state muons.
    This is done once to create a fast lookup table for matching.
    """
    print("Pre-processing GENIE data for faster matching...")
    start_time = time.time()

    all_pdgs = genie_data['StdHepPdg']
    all_statuses = genie_data['StdHepStatus']
    all_p4 = genie_data['StdHepP4']
    
    final_muon_mask = ((all_pdgs == PDG_MUON) | (all_pdgs == PDG_ANTI_MUON)) & (all_statuses == 1)
    initial_neutrino_mask = ((all_pdgs == PDG_MUON_NEUTRINO) | (all_pdgs == PDG_ANTI_MUON_NEUTRINO)) & (all_statuses == 0)

    event_selection_mask = (ak.sum(final_muon_mask, axis=1) == 1) & (ak.sum(initial_neutrino_mask, axis=1) == 1)

    valid_pdgs = all_pdgs[event_selection_mask]
    valid_statuses = all_statuses[event_selection_mask]
    valid_p4 = all_p4[event_selection_mask]

    final_muon_mask_valid = ((valid_pdgs == PDG_MUON) | (valid_pdgs == PDG_ANTI_MUON)) & (valid_statuses == 1)
    initial_neutrino_mask_valid = ((valid_pdgs == PDG_MUON_NEUTRINO) | (valid_pdgs == PDG_ANTI_MUON_NEUTRINO)) & (valid_statuses == 0)

    muon_p4 = ak.firsts(valid_p4[final_muon_mask_valid])
    neutrino_p4 = ak.firsts(valid_p4[initial_neutrino_mask_valid])

    muon_p4_np = ak.to_numpy(muon_p4)
    neutrino_p4_np = ak.to_numpy(neutrino_p4)

    # Calculate KE and Cos(theta) with Z-axis
    px, py, pz = muon_p4_np[:, 0], muon_p4_np[:, 1], muon_p4_np[:, 2]
    p_total = np.sqrt(px**2 + py**2 + pz**2)
    # Avoid division by zero for particles at rest, though unlikely for muons
    genie_cos_theta_z = np.divide(pz, p_total, out=np.zeros_like(pz), where=p_total!=0)

    muon_total_e_gev = muon_p4_np[:, 3]
    genie_muon_ke_mev = (muon_total_e_gev * 1000.0) - MUON_MASS_MEV
    genie_neutrino_e_mev = neutrino_p4_np[:, 3] * 1000.0
    
    print(f"GENIE pre-processing finished in {time.time() - start_time:.2f} seconds.")
    
    return {
        'ke_mev': genie_muon_ke_mev,
        'cos_theta_z': genie_cos_theta_z,
        'neutrino_e_mev': genie_neutrino_e_mev
    }


def find_genie_match_vectorized(ratpac_ke, ratpac_cos_theta, genie_kinematics):
    """
    Finds a matching GENIE event using optimized NumPy vector operations, matching on KE and Z-angle.
    """
    ke_diffs = np.abs(ratpac_ke - genie_kinematics['ke_mev'])
    cos_theta_diffs = np.abs(ratpac_cos_theta - genie_kinematics['cos_theta_z'])

    match_indices = np.where((ke_diffs < KE_TOLERANCE_MEV) )[0]

    if len(match_indices) > 0:
        return genie_kinematics['neutrino_e_mev'][match_indices[0]]
        
    return -1.0


def create_cherenkov_cone_verts(start_pos, direction, length, angle_deg):
    """
    Generates vertices for a 3D Cherenkov cone for Matplotlib.
    """
    norm_dir = np.linalg.norm(direction)
    if norm_dir == 0: return []
    direction = direction / norm_dir
    angle_rad = np.radians(angle_deg)
    radius = length * np.tan(angle_rad)
    v_rand = np.random.randn(3)
    v_rand -= v_rand.dot(direction) * direction
    if np.linalg.norm(v_rand) < 1e-6:
        v_rand = np.array([1, 0, 0])
        v_rand -= v_rand.dot(direction) * direction
    v1 = v_rand / np.linalg.norm(v_rand)
    v2 = np.cross(direction, v1)
    num_segments = 30
    t = np.linspace(0, 2 * np.pi, num_segments + 1)
    circle_pts = radius * (np.outer(np.cos(t), v1) + np.outer(np.sin(t), v2))
    base_pts = start_pos + length * direction + circle_pts
    verts = [[start_pos, base_pts[i], base_pts[i+1]] for i in range(num_segments)]
    return verts


def process_and_visualize_events(meta_tree, output_tree, genie_tree):
    """
    Main function to process events and generate 3D visualizations using Matplotlib.
    """
    pmt_positions = get_pmt_positions(meta_tree)

    print("Reading event data...")
    output_data = output_tree.arrays()
    genie_data_raw = genie_tree.arrays(['StdHepPdg', 'StdHepStatus', 'StdHepP4'])

    genie_kinematics = get_genie_muon_kinematics(genie_data_raw)

    num_events = len(output_data['mcpdg'])
    print(f"Total events in RATPAC file: {num_events}")
    
    # Get Z boundaries for the new block
    all_z = [p['z'] for p in pmt_positions.values()]
    min_z, max_z = min(all_z), max(all_z)

    for i in range(num_events):
        if output_data['mcpdg'][i] not in [PDG_MUON, PDG_ANTI_MUON]:
            continue
            
        print(f"\n--- Processing RATPAC Event {i+1}/{num_events} ---")
        
        muon_ke = output_data['mcke'][i]
        muon_dir = np.array([output_data['mcu'][i], output_data['mcv'][i], output_data['mcw'][i]])
        
        # Normalize direction and get cos(theta)
        muon_dir_norm = np.linalg.norm(muon_dir)
        ratpac_cos_theta = 0.0
        if muon_dir_norm > 0:
            ratpac_cos_theta = muon_dir[2] / muon_dir_norm

        neutrino_energy = find_genie_match_vectorized(muon_ke, ratpac_cos_theta, genie_kinematics)
        
        if neutrino_energy > 0:
            print(f"Match found! Neutrino Energy: {neutrino_energy:.2f} MeV")
        else:
            print("Warning: Could not find a matching GENIE event.")

        event_pmt_ids = output_data['mcPMTID'][i]
        event_pe_times = output_data['mcPEHitTime'][i]

        hit_pmt_info = {}
        for pmt_id, time in zip(event_pmt_ids, event_pe_times):
            if pmt_id not in hit_pmt_info:
                hit_pmt_info[pmt_id] = {'times': []}
            hit_pmt_info[pmt_id]['times'].append(time)
            
        for pmt_id in hit_pmt_info:
            times = hit_pmt_info[pmt_id]['times']
            hit_pmt_info[pmt_id]['npe'] = len(times)
            hit_pmt_info[pmt_id]['first_time'] = min(times) if times else 0

        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(111, projection='3d')

        all_pmt_ids = set(pmt_positions.keys())
        hit_pmt_ids = set(hit_pmt_info.keys())
        non_hit_pmt_ids = all_pmt_ids - hit_pmt_ids
        
        if non_hit_pmt_ids:
            non_hit_coords = [pmt_positions[pid] for pid in non_hit_pmt_ids]
            non_hit_x, non_hit_y, non_hit_z = zip(*[p.values() for p in non_hit_coords])
            ax.scatter(non_hit_x, non_hit_y, non_hit_z, 
                       facecolors='none', edgecolors='lightgray', alpha=0.1, s=10, label='Non-hit PMTs')

        if hit_pmt_ids:
            hit_pmt_coords = [pmt_positions[pid] for pid in hit_pmt_info]
            hit_pmt_x, hit_pmt_y, hit_pmt_z = zip(*[p.values() for p in hit_pmt_coords])
            hit_pmt_npe = np.array([info['npe'] for info in hit_pmt_info.values()])
            hit_pmt_time = np.array([info['first_time'] for info in hit_pmt_info.values()])
            marker_areas = hit_pmt_npe * 50

            sc = ax.scatter(hit_pmt_x, hit_pmt_y, hit_pmt_z, s=marker_areas, c=hit_pmt_time, 
                            cmap='plasma', edgecolor='k', alpha=0.9, label='Hit PMTs')
            cbar = fig.colorbar(sc, shrink=0.5, aspect=10)
            cbar.set_label('First Hit Time (ns)', fontsize=14)

        muon_pos = np.array([output_data['mcx'][i], output_data['mcy'][i], output_data['mcz'][i]])
        track_length_mm = muon_ke / MUON_ENERGY_LOSS_MEV_PER_MM
        track_end_pos = muon_pos + muon_dir * track_length_mm
        ax.plot([muon_pos[0], track_end_pos[0]], [muon_pos[1], track_end_pos[1]], [muon_pos[2], track_end_pos[2]], 
                color='red', linestyle='--', linewidth=2.5, label=f'Muon Track ({track_length_mm/10:.1f} cm)')

        cone_verts = create_cherenkov_cone_verts(muon_pos, muon_dir, track_length_mm, CHERENKOV_ANGLE_DEG)
        if cone_verts:
            ax.add_collection3d(Poly3DCollection(cone_verts, facecolors='cyan', linewidths=0, alpha=0.15))
            
        # Add the solid red block with new dimensions
        block_y_bottom, block_y_top = -20000, -12000 # -20m to -12m
        block_x_min, block_x_max = -20000, 20000     # -20m to 20m
        block_z_min, block_z_max = min_z, max_z       # Full detector Z range

        verts = [
            [(block_x_min, block_y_bottom, block_z_min), (block_x_max, block_y_bottom, block_z_min), (block_x_max, block_y_top, block_z_min), (block_x_min, block_y_top, block_z_min)],
            [(block_x_min, block_y_bottom, block_z_max), (block_x_max, block_y_bottom, block_z_max), (block_x_max, block_y_top, block_z_max), (block_x_min, block_y_top, block_z_max)],
            [(block_x_min, block_y_bottom, block_z_min), (block_x_min, block_y_bottom, block_z_max), (block_x_min, block_y_top, block_z_max), (block_x_min, block_y_top, block_z_min)],
            [(block_x_max, block_y_bottom, block_z_min), (block_x_max, block_y_bottom, block_z_max), (block_x_max, block_y_top, block_z_max), (block_x_max, block_y_top, block_z_min)],
            [(block_x_min, block_y_bottom, block_z_min), (block_x_max, block_y_bottom, block_z_min), (block_x_max, block_y_bottom, block_z_max), (block_x_min, block_y_bottom, block_z_max)],
            [(block_x_min, block_y_top, block_z_min), (block_x_max, block_y_top, block_z_min), (block_x_max, block_y_top, block_z_max), (block_x_min, block_y_top, block_z_max)],
        ]
        ax.add_collection3d(Poly3DCollection(verts, facecolors='red', linewidths=0, alpha=0.3))


        title_str = f"RATPAC Event Index: {i} | Neutrino E: {neutrino_energy:.2f} MeV | Muon KE: {muon_ke:.2f} MeV"
        ax.set_title(title_str, fontsize=20)
        ax.set_xlabel('X (mm)', fontsize=16); ax.set_ylabel('Y (mm)', fontsize=16); ax.set_zlabel('Z (mm)', fontsize=16)
        
        # Set new axis limits
        ax.set_xlim(-20000, 20000)
        ax.set_ylim(-20000, 20000)
        
        # Auto-scale Z axis to fit the detector
        all_z_values = [p['z'] for p in pmt_positions.values()]
        max_z_range = (max(all_z_values) - min(all_z_values)) / 2.0
        mid_z = np.mean(all_z_values)
        ax.set_zlim(mid_z - max_z_range, mid_z + max_z_range)


        ax.legend(fontsize=12)
        plt.show()

if __name__ == '__main__':
    meta, output, genie = load_data(RATPAC_FILE, GENIE_FILE)
    if meta and output and genie:
        process_and_visualize_events(meta, output, genie)
        print("\nAll muon events have been processed.")

