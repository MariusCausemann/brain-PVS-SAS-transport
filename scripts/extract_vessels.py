import nibabel
import pyvista as pv
import numpy as np
from skimage.filters import frangi, sato
import skimage.morphology as skim
from skimage.measure import label
import kimimaro
import os
import networkx as nx


def np2pv(arr, resolution, origin):
    grid = pv.ImageData(dimensions=arr.shape + np.array((1, 1, 1)),
                        spacing=resolution, origin=origin)
    grid[f"data"] = arr.flatten(order="F")
    return grid

def skel_to_mesh(skel):
    n_segs = skel.edges.shape[0]
    cells = np.concatenate([np.ones((n_segs,1), dtype=int)*2, skel.edges], axis=1)
    celltypes = np.full(n_segs, fill_value=pv.CellType.LINE, dtype=np.uint8)
    mesh = pv.UnstructuredGrid(cells.ravel(), celltypes, skel.vertices)
    mesh["radius"] = skel.radius
    return mesh

def skeletonize(img):
    skels = kimimaro.skeletonize(img, anisotropy=(1, 1, 1),
                             teasar_params={"scale": scale, "const": const,})
    joined_skels = kimimaro.join_close_components(skels.values(), radius=50)
    ds_skel = joined_skels.downsample(2)
    return ds_skel 

def as_networkx(skel, nroots):
    G = nx.DiGraph()
    G.add_edges_from(skel.edges)
    nx.set_node_attributes(G, {i:pos for i, pos in enumerate(skel.vertices)}, "pos")
    nx.set_node_attributes(G, {i:r for i, r in enumerate(skel.radius)}, "radius")
    degrees = np.array([G.degree(i) for i in range(G.number_of_nodes())])
    roots = np.where(degrees==1)[0]
    vec_dict = {}
    zcoords = skel.vertices[:,2]
    rs = roots[np.argsort(zcoords[roots])][:nroots]
    for e in nx.edge_bfs(G, source=rs, orientation="ignore"):
        vec = G.nodes[e[1]]["pos"] - G.nodes[e[0]]["pos"]
        orientation = 1 if e[2] == "forward" else -1
        vec_dict[e[:2]] = vec * orientation / np.linalg.norm(vec)
    nx.set_edge_attributes(G, vec_dict, "vector")
    nx.set_node_attributes(G, {i:1 if i in rs else 0 for i in range(G.number_of_nodes())}, "root")
    return G

def nx_to_pv(G):
    n_segs = G.number_of_edges()
    egdes = np.array(G.edges)
    points = np.array([G.nodes[i]["pos"] for i in range(G.number_of_nodes())])
    cells = np.concatenate([np.ones((n_segs,1), dtype=int)*2, egdes], axis=1)
    celltypes = np.full(n_segs, fill_value=pv.CellType.LINE, dtype=np.uint8)
    mesh = pv.UnstructuredGrid(cells.ravel(), celltypes, points)
    radii = np.array([G.nodes[i]["radius"] for i in range(G.number_of_nodes())])
    vec = np.array([data["vector"] for _, _, data in G.edges(data=True)])
    root = np.array([G.nodes[i]["root"] for i in range(G.number_of_nodes())])
    mesh["radius"] = radii
    mesh["orientation"] = vec
    mesh["root"] = root
    return mesh

scale = 1.0
const = 1

os.makedirs("../mesh/networks", exist_ok=True)

data = nibabel.load("../data/pcbi.1007073.s007.nii.gz")
img = data.get_fdata().astype(int)
arterial_skel = skeletonize(img)
G_art = as_networkx(arterial_skel, nroots=5)
arterial_network = nx_to_pv(G_art)
arterial_network.save("../mesh/networks/arteries.vtk")

data = nibabel.load("../data/pcbi.1007073.s008.nii.gz")
img = data.get_fdata().astype(int)
venous_skel = skeletonize(img)
G_ven = as_networkx(venous_skel, nroots=50)
venous_network = nx_to_pv(G_ven)
venous_network.save("../mesh/networks/venes.vtk")


