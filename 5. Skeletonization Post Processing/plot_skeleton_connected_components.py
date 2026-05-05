#Author: Lucas Philipp
#Subdivide skeleton into thick tubes and thin bridges by thresholding the max inscribed ball radii. 
#Uses connected components to further subdivide skeleton at junction nodes with connectivity degree >= 3.

import numpy as np
import networkx as nx
import napari
from matplotlib import cm

def load_off_mesh(path):
    with open(path, "r") as f:
        first_line = f.readline().strip()
        if first_line != "OFF":
            raise ValueError("Not a valid OFF file")
        counts = f.readline().split()
        n_vertices = int(counts[0])
        n_faces = int(counts[1])
        vertices = []
        for _ in range(n_vertices):
            x, y, z = map(float, f.readline().split())
            vertices.append([z, y, x])  # napari uses (z,y,x)
        vertices = np.array(vertices)
        faces = []
        for _ in range(n_faces):
            vals = list(map(int, f.readline().split()))
            n = vals[0]
            if n == 3:
                faces.append(vals[1:4])
        faces = np.array(faces)
    return vertices, faces

def split_at_junctions(G, min_vertices=10):
    """
    Remove all nodes with degree >= 3 from graph G, then return
    connected components of the pruned graph as a list of sets.
    Only components with >= min_vertices nodes are kept.
    """
    junction_nodes = {n for n, deg in G.degree() if deg >= 3}
    G_split = G.copy()
    G_split.remove_nodes_from(junction_nodes)
    return [
        comp for comp in nx.connected_components(G_split)
        if len(comp) >= min_vertices
    ]

def branch_to_lines(branch_vertices):
    lines = []
    for (v1, v2) in edges:
        if v1 in branch_vertices and v2 in branch_vertices:
            lines.append([vertices[v1], vertices[v2]])
    return np.array(lines) if lines else np.empty((0, 2, 3))


def longest_path_in_branch(branch_vertices, graph):
    """
    Compute contour length of branch.
    Returns (path_node_list, length_in_um)
    """
    subgraph = graph.subgraph(branch_vertices).copy()

    for u, v in subgraph.edges():
        dist = np.linalg.norm(vertices[u] - vertices[v])
        subgraph[u][v]['weight'] = dist

    start = next(iter(branch_vertices))

    # First Dijkstra pass: find farthest node from arbitrary start
    lengths1 = nx.single_source_dijkstra_path_length(subgraph, start, weight='weight')
    far_node1 = max(lengths1, key=lengths1.get)

    # Second Dijkstra pass: find farthest node from far_node1
    lengths2, paths2 = nx.single_source_dijkstra(subgraph, far_node1, weight='weight')
    far_node2 = max(lengths2, key=lengths2.get)

    longest_path = paths2[far_node2]
    longest_length = lengths2[far_node2]
    return longest_path, longest_length

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

# Load data

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

voxel_size = 4/1000 # [um]

# Scale skeleton
vertices *= voxel_size
radii *= voxel_size
#edges are edge indices

line_segments *= voxel_size
max_inscribed_ball_radii *= voxel_size

mesh_vertices *= voxel_size
#mesh_faces are indices

###

#kawagutii_run3_se2_right_cell
#threshold = 35*4/1000 #[um]
#min_vertices = 10

#kawagutii_run2_se2_top_right_cell
threshold = 30*4/1000 #[um]
min_vertices = 20

#kawagutii_run2_se2_middle_cell
#threshold = 33*4/1000 #[um]
#min_vertices = 30

###

# Build graph
G_low = nx.Graph() #radii below threshold
G_high = nx.Graph() #radii above threshold

#run plot_skeleton.py first
for (v1, v2), r in zip(edges, max_inscribed_ball_radii):
    if r < threshold:
        G_low.add_edge(v1, v2)
    else:
        G_high.add_edge(v1, v2)


low_branches = [
    comp for comp in nx.connected_components(G_low)
    if len(comp) >= min_vertices
]

high_branches = [
    comp for comp in nx.connected_components(G_high)
    if len(comp) >= min_vertices
]

# Split skeleton at junction nodes (degree >= 3)
###
#uncomment indices_to_remove to prune connected components based on visual inspection. these indices only apply for the threshold & min_vertices values indicated
###

low_branches = split_at_junctions(G_low, min_vertices=min_vertices)
#kawagutii_run3_se2_right_cell
#threshold = 35*4/1000 [um]
#min_vertices = 10
# indices_to_remove = {
#     6, 8, 10, 11, 15, 17, 21, 22, 23, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,
#     38, 39, 41, 42, 43, 45, 46, 48, 49, 50, 51, 53, 54, 55, 56, 57, 59, 60, 61, 63,
#     64, 66, 69, 70, 73, 74, 76, 77, 78, 80, 81, 83, 84, 86, 87, 90, 91, 92, 93, 94,
#     96, 97, 98, 99, 101, 104, 106, 108, 109, 116, 117, 118, 120, 122, 123, 125, 130,
#     138, 143, 145, 146, 147, 149, 153, 154, 155, 156, 165, 166, 170, 171, 173, 174,
#     176, 178, 179, 180, 185, 186, 187, 188, 189, 193, 194, 195, 196, 197, 199, 200,
#     201, 205, 206, 207, 209, 210, 211, 212, 217, 218, 219, 220, 224, 225, 228, 229,
#     231, 233, 234, 235, 237, 242, 244, 249, 250, 252, 253, 254, 255, 258, 259, 260,
#     261, 262, 263, 265, 267, 268, 269, 271, 275, 276, 280, 281, 282, 283, 285, 286,
#     287, 288, 289, 290, 291, 292, 293, 294
# }

#kawagutii_run2_se2_top_right_cell
#threshold = 30*4/1000 [um]
#min_vertices = 20
indices_to_remove = {0, 1, 2, 3, 4, 11, 12, 13, 16, 19, 22, 23, 24, 27, 29, 32, 37,
38, 39, 40, 41, 43, 44, 46, 48, 49, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 65, 66, 
68, 70, 71, 72, 74, 75, 81, 82, 83, 84, 85, 86, 88, 91, 93, 121, 122, 127, 132, 134, 136, 
138, 139, 140, 148, 150, 151, 160, 161, 162, 163, 165, 167, 175, 176, 177, 179, 180, 181, 
188, 189, 197, 200, 205, 206, 207, 208, 212, 225, 228, 233, 244, 246, 248, 249, 250, 251, 252
}

#kawagutii_run2_se2_middle_cell
#threshold = 33*4/1000 [um]
#min_vertices = 30
# indices_to_remove = {0, 1, 2, 4, 7, 10, 12, 13, 14, 16, 17, 20, 23, 26, 31, 32, 33, 40, 41, 42, 43, 
# 47, 52, 54, 56, 57, 58, 61, 63, 64, 65, 67, 71, 73, 74, 75, 78, 86, 87, 91, 94, 
# 95, 98, 99, 108, 111, 114, 117, 118, 122, 124, 125, 126, 127, 128, 131, 144, 148, 
# 149, 145, 146, 147, 150, 153, 158, 163, 164, 165, 176, 181, 182, 185, 186, 188, 
# 192, 199, 208, 213, 214, 215, 216, 217, 220, 221, 222, 223, 224, 225, 226, 227, 
# 228, 230
# }

low_branches = [b for i, b in enumerate(low_branches) if i not in indices_to_remove]

#if you want to visualize just a subset of the branches (loads faster)
low_branches = low_branches[:10]
#low_branches = low_branches[10:20]
#low_branches = low_branches[20:30]
#low_branches = low_branches[30:40]
#low_branches = low_branches[40:50]
#low_branches = low_branches[50:60]
#low_branches = low_branches[60:70]
#low_branches = low_branches[70:80]
#low_branches = low_branches[80:90]
#low_branches = low_branches[90:100]
#low_branches = low_branches[100:110]
#low_branches = low_branches[110:120]
#low_branches = low_branches[120:130]
#low_branches = low_branches[130:140]
#low_branches = low_branches[140:150]
#low_branches = low_branches[150:160]
#low_branches = low_branches[160:170]
#low_branches = low_branches[170:180]
#low_branches = low_branches[180:190]
#low_branches = low_branches[190:200]
#low_branches = low_branches[200:210]
#low_branches = low_branches[210:220]
#low_branches = low_branches[220:230]
#low_branches = low_branches[230:240]
#low_branches = low_branches[240:250]
#low_branches = low_branches[250:260]
#low_branches = low_branches[260:270]
#low_branches = low_branches[270:280]
#low_branches = low_branches[280:290]

print("Number of branches below radii threshold:", len(low_branches))
print("Number of branches above radii threshold:", len(high_branches))

# Compute longest paths for all high or low branches

print("\nComputing longest paths for branches below radii threshold...")
branch_longest_paths_low = []
branch_longest_paths_high = []

for i, branch in enumerate(low_branches):
    path, length = longest_path_in_branch(branch, G_low)
    branch_longest_paths_low.append({
        'branch_idx': i,
        'branch_vertices': branch,
        'path': path,
        'length_um': length,
        'n_vertices': len(branch)
    })
    print(f"  Branch {i:3d} | vertices: {len(branch):5d} | "
          f"longest path: {length:.4f} um")

print("\nComputing longest paths for branches above radii threshold...")
for i, branch in enumerate(high_branches):
    path, length = longest_path_in_branch(branch, G_high)
    branch_longest_paths_high.append({
        'branch_idx': i,
        'branch_vertices': branch,
        'path': path,
        'length_um': length,
        'n_vertices': len(branch)
    })
    print(f"  Branch {i:3d} | vertices: {len(branch):5d} | "
          f"longest path: {length:.4f} um | hops: {len(path) - 1}")

viewer = napari.Viewer(ndisplay=3)

viewer.add_surface(
    (mesh_vertices, mesh_faces),
    colormap="gray",
    opacity=0.3,
    name="Volume Surface",
    blending="translucent_no_depth",
)

# Low branches (shades of blue)
blues = cm.Blues(np.linspace(0.3, 1, len(low_branches)))
for branch, color in zip(low_branches, blues):
    lines = branch_to_lines(branch)
    if len(lines) > 0:
        viewer.add_shapes(
            lines,
            shape_type='line',
            edge_color=[color],
            edge_width=5/1000,
            name='low_branch'
        )

# # High branches (shades of red)
# reds = cm.Reds(np.linspace(0.3, 1, len(high_branches)))
# for branch, color in zip(high_branches, reds):
#     lines = branch_to_lines(branch)
#     if len(lines) > 0:
#         viewer.add_shapes(
#             lines,
#             shape_type='line',
#             edge_color=[color],
#             edge_width=5 / 1000,
#             name='high_branches'
#         )

viewer.theme = "light"

viewer.scale_bar.visible = True
viewer.scale_bar.unit = "µm" # or "nm"
viewer.scale_bar.position = "bottom_right"
viewer.scale_bar.color = "black"  # better for white background
viewer.scale_bar.length = 1

napari.run()

viewer.camera.zoom *= 1.4   # 40% zoom in
napari.run()