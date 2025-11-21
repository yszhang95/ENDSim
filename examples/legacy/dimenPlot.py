import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import argparse
import matplotlib

# --- Configuration & Calculations ---

# Default space size
DEFAULT_SPACE_SIZE = 100.0 # meters

# Rectangle areas
SMALL_RECT_AREA = 1.5 # m^2
LARGE_RECT_AREA = 4.5 # m^2
NUM_SMALL_RECTS = 4
NUM_LARGE_RECTS = 18

# Colors (matching the JS version as closely as possible)
COLORS = {
    'xy': {'large': '#0000ff', 'small': '#add8e6'}, # Blue / Light Blue
    'yz': {'large': '#ff0000', 'small': '#ffc0cb'}, # Red / Pink
    'xz': {'large': '#008000', 'small': '#90ee90'}  # Green / Light Green
}
BOX_COLOR = '#333333' # Bounding box color

# Gap between tiles
GAP = 0.05

def calculate_dimensions(plane_type):
    """Calculates rectangle dimensions based on plane type and rules."""
    if plane_type == 'xy': # Long leg Y (Lx <= Ly)
        small_dims = {'width': 1.0, 'height': 1.5} # width=X, height=Y
        large_dims = {'width': 1.5, 'height': 3.0} # width=X, height=Y
    elif plane_type == 'yz': # Long leg Z (Ly <= Lz)
        small_dims = {'width': 1.0, 'height': 1.5} # width=Y, height=Z
        large_dims = {'width': 1.5, 'height': 3.0} # width=Y, height=Z
    elif plane_type == 'xz': # Long leg Z (Lx <= Lz)
        small_dims = {'width': 1.0, 'height': 1.5} # width=X, height=Z
        large_dims = {'width': 1.5, 'height': 3.0} # width=X, height=Z
    else:
        raise ValueError("Invalid plane_type. Must be 'xy', 'yz', or 'xz'.")
    return small_dims, large_dims

def create_rectangle_vertices(center_x, center_y, center_z, width, height, plane_type):
    """Creates the 4 vertices for a rectangle on a specific plane."""
    half_w = width / 2.0
    half_h = height / 2.0

    if plane_type == 'xy': # Plane at Z = center_z
        # Vertices order: bottom-left, bottom-right, top-right, top-left (counter-clockwise)
        verts = [
            (center_x - half_w, center_y - half_h, center_z),
            (center_x + half_w, center_y - half_h, center_z),
            (center_x + half_w, center_y + half_h, center_z),
            (center_x - half_w, center_y + half_h, center_z)
        ]
    elif plane_type == 'yz': # Plane at X = center_x
        # Width corresponds to Y, Height corresponds to Z
        verts = [
            (center_x, center_y - half_w, center_z - half_h),
            (center_x, center_y + half_w, center_z - half_h),
            (center_x, center_y + half_w, center_z + half_h),
            (center_x, center_y - half_w, center_z + half_h)
        ]
    elif plane_type == 'xz': # Plane at Y = center_y
        # Width corresponds to X, Height corresponds to Z
        verts = [
            (center_x - half_w, center_y, center_z - half_h),
            (center_x + half_w, center_y, center_z - half_h),
            (center_x + half_w, center_y, center_z + half_h),
            (center_x - half_w, center_y, center_z + half_h)
        ]
    else:
         raise ValueError("Invalid plane_type.")

    return [verts] # Return as a list containing one polygon

def create_tiled_plane_polygons(plane_type):
    """Generates a list of polygons (vertices) for the tiled plane."""
    small_rect_dims, large_rect_dims = calculate_dimensions(plane_type)
    plane_colors = COLORS[plane_type]
    polygons = []
    colors = []

    # --- Define layout strategy (same as JS version) ---
    large_grid_width_dim = 3 * (large_rect_dims['width'] + GAP) - GAP
    large_grid_height_dim = 6 * (large_rect_dims['height'] + GAP) - GAP
    small_grid_width_dim = 2 * (small_rect_dims['width'] + GAP) - GAP
    small_grid_height_dim = 2 * (small_rect_dims['height'] + GAP) - GAP

    total_width_dim = large_grid_width_dim + small_grid_width_dim + GAP
    total_height_dim = max(large_grid_height_dim, small_grid_height_dim)

    start_w = -total_width_dim / 2.0  # Start coordinate for the 'width' dimension on the plane
    start_h = -total_height_dim / 2.0 # Start coordinate for the 'height' dimension on the plane

    # Create Large Rectangles (18) - 3 columns, 6 rows
    current_h = start_h
    for row in range(6):
        current_w = start_w
        for col in range(3):
            center_w = current_w + large_rect_dims['width'] / 2.0
            center_h = current_h + large_rect_dims['height'] / 2.0

            # Map center_w, center_h to x, y, z based on plane_type
            if plane_type == 'xy':
                verts = create_rectangle_vertices(center_w, center_h, 0, large_rect_dims['width'], large_rect_dims['height'], plane_type)
            elif plane_type == 'yz': # width=Y, height=Z
                verts = create_rectangle_vertices(0, center_w, center_h, large_rect_dims['width'], large_rect_dims['height'], plane_type)
            elif plane_type == 'xz': # width=X, height=Z
                verts = create_rectangle_vertices(center_w, 0, center_h, large_rect_dims['width'], large_rect_dims['height'], plane_type)

            polygons.extend(verts)
            colors.append(plane_colors['large'])
            current_w += large_rect_dims['width'] + GAP
        current_h += large_rect_dims['height'] + GAP

    # Create Small Rectangles (4) - 2 columns, 2 rows, next to the large grid
    current_w_small_start = start_w + large_grid_width_dim + GAP
    current_h = start_h
    for row in range(2):
        current_w = current_w_small_start
        for col in range(2):
            center_w = current_w + small_rect_dims['width'] / 2.0
            center_h = current_h + small_rect_dims['height'] / 2.0

            # Map center_w, center_h to x, y, z based on plane_type
            if plane_type == 'xy':
                verts = create_rectangle_vertices(center_w, center_h, 0, small_rect_dims['width'], small_rect_dims['height'], plane_type)
            elif plane_type == 'yz': # width=Y, height=Z
                verts = create_rectangle_vertices(0, center_w, center_h, small_rect_dims['width'], small_rect_dims['height'], plane_type)
            elif plane_type == 'xz': # width=X, height=Z
                verts = create_rectangle_vertices(center_w, 0, center_h, small_rect_dims['width'], small_rect_dims['height'], plane_type)

            polygons.extend(verts)
            colors.append(plane_colors['small'])
            current_w += small_rect_dims['width'] + GAP
        current_h += small_rect_dims['height'] + GAP

    return polygons, colors

def draw_bounding_box(ax, size):
    """Draws the wireframe bounding box."""
    half_size = size / 2.0
    points = np.array([
        [-half_size, -half_size, -half_size],
        [ half_size, -half_size, -half_size],
        [ half_size,  half_size, -half_size],
        [-half_size,  half_size, -half_size],
        [-half_size, -half_size,  half_size],
        [ half_size, -half_size,  half_size],
        [ half_size,  half_size,  half_size],
        [-half_size,  half_size,  half_size]
    ])
    edges = [
        [points[0], points[1]], [points[1], points[2]], [points[2], points[3]], [points[3], points[0]], # bottom face
        [points[4], points[5]], [points[5], points[6]], [points[6], points[7]], [points[7], points[4]], # top face
        [points[0], points[4]], [points[1], points[5]], [points[2], points[6]], [points[3], points[7]]  # side edges
    ]
    for edge in edges:
        line = np.array(edge)
        ax.plot(line[:, 0], line[:, 1], line[:, 2], color=BOX_COLOR, alpha=0.6)


# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize a 3D space with a tiled plane.")
    parser.add_argument(
        "--size",
        type=float,
        default=DEFAULT_SPACE_SIZE,
        help=f"Side length of the cubic space (default: {DEFAULT_SPACE_SIZE}m)."
    )
    parser.add_argument(
        "--plane",
        type=str,
        choices=['xy', 'yz', 'xz'],
        default='xy',
        help="Which plane to display (xy, yz, or xz) (default: xy)."
    )
    args = parser.parse_args()

    space_size = args.size
    plane_to_show = args.plane

    # --- Plotting Setup ---
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Generate polygons for the selected plane
    polygons, colors = create_tiled_plane_polygons(plane_to_show)

    # Create the 3D polygon collection
    poly_collection = Poly3DCollection(polygons, facecolors=colors, edgecolors='k', linewidths=0.5, alpha=0.9)
    ax.add_collection3d(poly_collection)

    # Draw the bounding box
    draw_bounding_box(ax, space_size)

    # --- Axes Configuration ---
    half_size = space_size / 2.0
    plot_limit = half_size * 1.1 # Add a small margin

    ax.set_xlim([-plot_limit, plot_limit])
    ax.set_ylim([-plot_limit, plot_limit])
    ax.set_zlim([-plot_limit, plot_limit])

    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")

    # Set aspect ratio to equal
    # ax.set_aspect('equal', adjustable='box') # This causes error in older Matplotlib versions
    try:
        # Recommended method for Matplotlib >= 3.3
        ax.set_box_aspect([1,1,1])
    except AttributeError:
        # Fallback for older versions (might not be perfectly equal)
        print("Warning: Matplotlib version might be older than 3.3. "
              "Using manual limits for aspect ratio approximation.")
        # Calculate the range of data
        all_verts = np.concatenate([p[0] for p in polygons]) # Get all vertices
        box_points = np.array([ # Include box corners for limits
             [-half_size, -half_size, -half_size], [ half_size,  half_size,  half_size]
        ])
        all_points = np.vstack([all_verts, box_points])
        max_range = np.array([all_points[:,0].max()-all_points[:,0].min(),
                              all_points[:,1].max()-all_points[:,1].min(),
                              all_points[:,2].max()-all_points[:,2].min()]).max() / 2.0

        mid_x = (all_points[:,0].max()+all_points[:,0].min()) * 0.5
        mid_y = (all_points[:,1].max()+all_points[:,1].min()) * 0.5
        mid_z = (all_points[:,2].max()+all_points[:,2].min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        print(f"Adjusted limits to approx. equal aspect: +/- {max_range:.2f} around center.")


    ax.set_title(f"3D Space ({space_size}m Cube) showing {plane_to_show.upper()} Plane")

    plt.tight_layout()
    plt.show()

