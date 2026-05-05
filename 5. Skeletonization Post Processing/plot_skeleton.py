#Author: Lucas Philipp
#load output from Voxel Cores/Erosion Thickness and visualize radii of maximally inscribed balls
#color skeleton according to radii threshold

import numpy as np
import napari
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

def load_ply_skeleton(path):
    with open(path, "r") as f:
        lines = f.readlines()

    # --- Parse header ---
    i = 0
    n_vertices = 0
    n_edges = 0

    while not lines[i].startswith("end_header"):
        if lines[i].startswith("element vertex"):
            n_vertices = int(lines[i].split()[-1])
        if lines[i].startswith("element edge"):
            n_edges = int(lines[i].split()[-1])
        i += 1

    header_end = i + 1

    # --- Parse vertices ---
    vertices_raw = lines[header_end:header_end + n_vertices]
    vertices = []
    radii = []

    for line in vertices_raw:
        vals = list(map(float, line.split()))
        # Format: bt2 radius x y z
        radius = vals[1]
        x, y, z = vals[2:5] #2, 3, 4

        vertices.append([z, y, x])  # napari expects (z,y,x)
        radii.append(radius)

    vertices = np.array(vertices)
    radii = np.array(radii)

    # --- Parse edges ---
    edge_start = header_end + n_vertices
    edges_raw = lines[edge_start:edge_start + n_edges]

    edges = []
    for line in edges_raw:
        v1, v2 = map(int, line.split())
        edges.append([v1, v2])

    edges = np.array(edges)

    return vertices, edges, radii

def build_edge_line_segments(vertices, edges, radii):
    line_segments = []
    max_inscribed_ball_radii = []

    for v1, v2 in edges:
        line_segments.append([vertices[v1], vertices[v2]])
        max_inscribed_ball_radii.append((radii[v1] + radii[v2]) / 2.0)

    return np.array(line_segments), np.array(max_inscribed_ball_radii)

def load_off_mesh(path):
    with open(path, "r") as f:
        first_line = f.readline().strip()
        if first_line != "OFF":
            raise ValueError("Not a valid OFF file")

        counts = f.readline().split()
        n_vertices = int(counts[0])
        n_faces = int(counts[1])

        # Load vertices
        vertices = []
        for _ in range(n_vertices):
            x, y, z = map(float, f.readline().split())
            vertices.append([z, y, x])  # napari uses (z,y,x)

        vertices = np.array(vertices)

        # Load faces
        faces = []
        for _ in range(n_faces):
            vals = list(map(int, f.readline().split()))
            n = vals[0]
            if n == 3:  # triangular mesh
                faces.append(vals[1:4])

        faces = np.array(faces)

    return vertices, faces

#Load data

#load segmentation
#off_path = "/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/skeleton/kawagutii_run2_se2_middle_cell.off"
off_path = "/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/skeleton/kawagutii_run2_se2_top_right_cell.off"
#off_path = "/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/skeleton/kawagutii_run3_se2_right_cell.off"
mesh_vertices, mesh_faces = load_off_mesh(off_path)

#load skeleton
#ply_path = "/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/skeleton/kawagutii_run2_se2_middle_cell_thinned0.006_skel.ply"
ply_path = "/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/skeleton/kawagutii_run2_se2_top_right_cell_thinned0.006_skel.ply"
#ply_path = "/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/skeleton/kawagutii_run3_se2_right_cell_thinned0.006_skel.ply"

vertices, edges, radii = load_ply_skeleton(ply_path)
line_segments, max_inscribed_ball_radii = build_edge_line_segments(vertices, edges, radii) #napari expects line line_segments. line line_segments are a pair of xyz coordinates. will color each line segment according to max_inscribed_ball_radii

voxel_size = 4/(1000) # [um]

# Scale OFF mesh
mesh_vertices *= voxel_size

# Scale skeleton
vertices *= voxel_size
radii *= voxel_size
#edges are edge indices

line_segments *= voxel_size
max_inscribed_ball_radii *= voxel_size

#norm = Normalize(vmin=max_inscribed_ball_radii.min(), vmax=max_inscribed_ball_radii.max())
norm = Normalize(vmin=0/1000, vmax=300/1000)
cmap = cm.get_cmap("viridis")
radii_colormap = cmap(norm(max_inscribed_ball_radii))

threshold = 33*4/1000 #[um]
#threshold = 30*4/1000 #[um]
#threshold = 35*4/1000 #[um]

thresholded_colormap = np.where(max_inscribed_ball_radii < threshold, "blue", "red")

#plot histogram of radii
plt.hist(max_inscribed_ball_radii, bins=30, density=True)
plt.xlabel("Radius, r, [um]")
plt.ylabel("Count")
plt.show()

#np.savetxt("kawagutii_run2_se2_middle_cell_radius.csv", max_inscribed_ball_radii, delimiter=",", fmt="%.6f")
#np.savetxt("kawagutii_run2_se2_top_right_cell_radius.csv", max_inscribed_ball_radii, delimiter=",", fmt="%.6f")
#np.savetxt("kawagutii_run3_se2_right_cell_radius.csv", max_inscribed_ball_radii, delimiter=",", fmt="%.6f")

#Napari viewer
viewer = napari.Viewer(ndisplay=3)

viewer.add_shapes(
    line_segments,
    shape_type="line",
    edge_color=radii_colormap,
    edge_width=10/1000,
    name="Max inscribed ball radii",
)

viewer.add_shapes(
    line_segments,
    shape_type="line",
    edge_color=thresholded_colormap,
    edge_width=10/1000,
    name="Thresholded radii",
)

viewer.add_surface(
    (mesh_vertices, mesh_faces),
    colormap="gray",
    opacity=0.3,
    name="Volume Surface",
)

viewer.theme = "light"

viewer.scale_bar.visible = True
viewer.scale_bar.unit = "µm" # or "nm"
viewer.scale_bar.position = "bottom_right"
viewer.scale_bar.color = "black"
viewer.scale_bar.length = 1

napari.run()

viewer.camera.zoom *= 1.4   # 40% zoom in

# Plot horizontal colorbar
###
# Convert to nm
max_inscribed_ball_radii_nm = max_inscribed_ball_radii * 1000  # µm → nm
norm_nm = Normalize(vmin=0, vmax=300)

sm = ScalarMappable(cmap=cmap, norm=norm_nm)
sm.set_array([])

fig, ax = plt.subplots(figsize=(6, 1))
fig.subplots_adjust(bottom=0.5)

cbar = fig.colorbar(sm, cax=ax, orientation='horizontal')
cbar.set_label("Radius [nm]", fontsize=24)

# Set ticks in nm
cbar.set_ticks(np.linspace(0, 300, 7))
cbar.ax.tick_params(labelsize=24)

plt.show()
###