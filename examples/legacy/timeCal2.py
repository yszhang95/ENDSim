import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
# Speed of light in m/ns (meters per nanosecond)
c_m_ns = 0.2

# --- Formula for this geometry ---
def calculate_delta_t_parallel_pmts(L_perpendicular, d_parallel, c_val):
    """
    Calculates delta_t for PMTs parallel to muon track.
    delta_t is in ns if d_parallel is in m and c_val is in m/ns.
    L_perpendicular (perpendicular distance, your "y distance") is not used
    in the formula for this specific geometry but is included as an argument
    to match the meshgrid inputs.
    """
    # For PMTs parallel to the track, delta_t = d_parallel / c_val
    # It does not depend on L_perpendicular.
    return d_parallel / c_val

# --- Grid Setup ---
# Define ranges for L (perpendicular "y distance") and d_z ("dz") in meters
L_min, L_max, L_points = 0.0, 10.0, 50  # Perpendicular distance "y distance"
d_z_min, d_z_max, d_z_points = 0.0, 10.0, 50 # Parallel separation "dz"

L_range = np.linspace(L_min, L_max, L_points)
d_z_range = np.linspace(d_z_min, d_z_max, d_z_points)

# Create a grid of L and d_z values
# L_grid will vary along columns, d_z_grid will vary along rows
L_grid, d_z_grid = np.meshgrid(L_range, d_z_range)

# --- Calculate Delta_t on Grid ---
delta_t_grid = calculate_delta_t_parallel_pmts(L_grid, d_z_grid, c_m_ns)

# --- Plotting ---
plt.figure(figsize=(9, 7.5))

# Create a filled contour plot
# Choose contour levels (e.g., 20 levels, or specify explicitly)
# Manually set levels if desired for better control
levels = np.linspace(delta_t_grid.min(), delta_t_grid.max(), 25)
if delta_t_grid.min() == delta_t_grid.max(): # Handle case where d_z range results in single delta_t
    if delta_t_grid.min() == 0:
        levels = np.linspace(0, 1, 10) # Provide a default if all is zero
    else:
        levels = np.linspace(delta_t_grid.min()*0.9, delta_t_grid.max()*1.1, 10)


contour = plt.contourf(L_grid, d_z_grid, delta_t_grid, levels=levels, cmap='viridis', extend='max')

# Add contour lines for specific values for clarity
# Show fewer lines if the range is small or they become too crowded
if len(contour.levels) > 3:
    contour_lines_levels = contour.levels[::3] # Show every 3rd level line
else:
    contour_lines_levels = contour.levels

contour_lines = plt.contour(L_grid, d_z_grid, delta_t_grid, levels=contour_lines_levels, colors='white', linewidths=0.7)
plt.clabel(contour_lines, inline=True, fontsize=9, fmt='%.1f ns')

# Add a color bar
cbar = plt.colorbar(contour)
cbar.set_label('Time Difference $\Delta t = d_z / c$ (ns)')

# Labels and Title
plt.xlabel('Perpendicular Distance of Muon Track, $L$ (m) ', fontsize=14)
plt.ylabel('Separation between PMTs along Muon Path, $d_z$ (m)', fontsize=14)
plt.title('Cherenkov First Light $\Delta t$ (PMTs Parallel to Muon Track)', fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, linestyle=':', alpha=0.6)

# Show the plot
plt.show()
