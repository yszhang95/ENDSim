import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
n_water = 1.33
# Speed of light in m/ns (meters per nanosecond)
c_m_ns = 0.2

# --- Formula ---
def calculate_delta_t(L, d, n, c):
    """
    Calculates delta_t in ns if L, d are in m, c is in m/ns.
    L = distance of track from PMT1
    d = distance between PMT1 and PMT2
    n = refractive index
    c = speed of light in m/ns
    """
    if n <= 1:
        # No Cherenkov light if n <= 1
        return np.zeros_like(L + d) # Return array of zeros with correct shape

    n_factor = np.sqrt(n**2 - 1)

    # Ensure inputs are numpy arrays for element-wise operations
    L = np.asarray(L)
    d = np.asarray(d)

    # Calculate the geometric term
    # Use np.hypot(L, d) for potentially better numerical stability than sqrt(L**2+d**2)
    geometric_term = np.hypot(L, d) - L

    delta_t = (n_factor / c) * geometric_term
    return delta_t

# --- Grid Setup ---
# Define ranges for L and d in meters
L_min, L_max, L_points = 0, 10, 100
d_min, d_max, d_points = 0, 10, 100

L_range = np.linspace(L_min, L_max, L_points)
d_range = np.linspace(d_min, d_max, d_points)

# Create a grid of L and d values
L_grid, d_grid = np.meshgrid(L_range, d_range)

# --- Calculate Delta_t on Grid ---
delta_t_grid = calculate_delta_t(L_grid, d_grid, n_water, c_m_ns)

# --- Plotting ---
#plt.style.use('seaborn-v0_8-colorblind') # Nicer default style
plt.figure(figsize=(9, 7.5))

# Create a filled contour plot
# Choose contour levels (e.g., 20 levels, or specify explicitly)
levels = np.linspace(delta_t_grid.min(), delta_t_grid.max(), 25)
contour = plt.contourf(L_grid, d_grid, delta_t_grid, levels=levels, cmap='plasma', extend='max')

# Add contour lines for specific values for clarity
contour_lines = plt.contour(L_grid, d_grid, delta_t_grid, levels=contour.levels[::3], colors='black', linewidths=0.6) # Show every 3rd level line
plt.clabel(contour_lines, inline=True, fontsize=9, fmt='%.1f ns')

# Add a color bar
cbar = plt.colorbar(contour)
cbar.set_label('Time Difference $\Delta t = t_2 - t_1$ (ns)')

# Labels and Title
plt.xlabel('Distance of Muon Track from PMT1, $L$ (m)', fontsize=14)
plt.ylabel('Distance between PMTs, $d$ (m)', fontsize=14)
plt.title(f'Cherenkov First Light $\Delta t$ in Water (n={n_water:.2f}) vs. Geometry', fontsize=14)

# Optionally, make axis tick labels larger
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Set aspect ratio of the plot axes to be equal since both are meters
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True, linestyle=':', alpha=0.6)

# Show the plot
plt.show()
