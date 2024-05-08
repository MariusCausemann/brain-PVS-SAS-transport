import nibabel
import pyvista as pv
import numpy as np
from skimage.filters import frangi, sato
import skimage.morphology as skim
from skimage.measure import label
import kimimaro
import os
import networkx as nx
from cloudvolume import Bbox


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

def skeletonize(img, resolution):
    skels = kimimaro.skeletonize(img, anisotropy=resolution,
                             teasar_params={"scale": scale, "const": const,})
    joined_skels = kimimaro.join_close_components(skels.values(), radius=50)
    ds_skel = joined_skels.downsample(2)
    return ds_skel 

def generate_splines(G, point_ratio=1):
    G = G.to_undirected()
    degrees = np.array([G.degree(i) for i in range(G.number_of_nodes())])
    bifurcations = np.where(degrees > 2)[0]
    points = np.array([G.nodes[i]["pos"] for i in range(G.number_of_nodes())])
    Gwb = G.copy()
    Gwb.remove_nodes_from(bifurcations)
    splines = []

    for seg in list(nx.connected_components(Gwb)):
        g = G.subgraph(seg)
        roots = [node for (node, order) in g.degree() if order<=1]
        seg_nodes = seg
        for r in roots:
            seg_nodes = seg_nodes.union(set(G.neighbors(r)))
        g = G.subgraph(seg_nodes)
        newroots = [node for (node, order) in g.degree() if order==1]
        spl = pv.Spline(points=points[list(nx.dfs_postorder_nodes(g, source=newroots[0]))],
                         n_points=int(len(seg)*point_ratio))
        splines.append(spl)
    
    for bi in bifurcations:
        for nbi in G.neighbors(bi):
            if nbi in bifurcations:
                spl = pv.Spline(points=points[[bi, nbi]], n_points=int(2*point_ratio))
                splines.append(spl)

    return splines

def as_networkx(skel, nroots):
    G = nx.DiGraph()
    G.add_edges_from(skel.edges)
    nx.set_node_attributes(G, {i:pos for i, pos in enumerate(skel.vertices)}, "pos")
    nx.set_node_attributes(G, {i:r for i, r in enumerate(skel.radius)}, "radius")
    degrees = np.array([G.degree(i) for i in range(G.number_of_nodes())])

    # roots are all nodes that only have one connection/edge (degree = 1)
    roots = np.where(degrees==1)[0]
    vec_dict = {}
    zcoords = skel.vertices[:,2]

    # rs are the nroots roots with the nroots lowest z-coordinate (so actual roots)
    rs = roots[np.argsort(zcoords[roots])][:nroots]
    for e in nx.edge_bfs(G, source=rs, orientation="ignore"):
        vec = G.nodes[e[1]]["pos"] - G.nodes[e[0]]["pos"]
        orientation = 1 if e[2] == "forward" else -1
        vec_dict[e[:2]] = vec * orientation / np.linalg.norm(vec)
    nx.set_edge_attributes(G, vec_dict, "vector")

    # MER: I tried to make sense of the above and below, please edit at will.
    ROOT = 2
    EXITS = 1
    NOTROOT = 0
    def root_mark(i):
        if i in rs: return ROOT
        if i in roots: return EXITS
        return NOTROOT

    nx.set_node_attributes(G, {i:root_mark(i)
                               for i in range(G.number_of_nodes())}, "root")
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

def plot_radii(splines, filename):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrows=int(len(splines)/3),
                             ncols=3, sharey=True,
                             sharex=True, figsize=(5, int(len(splines)/4)))
    
    for ax,spl in zip(axes.flatten(),splines):
        radii = spl["radius"]
        arcl = spl["arc_length"]
        ax.plot(arcl, radii, "-*")

    fig.supxlabel("vessel segment length (mm)", y=0.001)
    fig.supylabel("radius (mm)")
    fig.suptitle("radius over segment length", y=0.999)
    plt.tight_layout()
    plt.savefig(filename)

def as_tubes(splines):
    tubes = pv.MultiBlock()
    for spl in splines:
        tubes.append(spl.tube(scalars="radius", absolute=True))
    return tubes


scale = 0.5
const = 0.5

os.makedirs("../mesh/networks", exist_ok=True)
os.makedirs("../plots", exist_ok=True)

files = ["../data/pcbi.1007073.s007.nii.gz", "../data/pcbi.1007073.s008.nii.gz"]
names = ["arteries", "venes"]
nroots = [3, 50]

for file, name, nr in zip(files, names, nroots):
    data = nibabel.load(file)
    resolution = data.header["pixdim"][1:4]
    img = data.get_fdata().astype(int)
    skel = skeletonize(img, resolution)
    #if name=="arteries":
    #    skel = skel.crop(Bbox([0, 0, 78], img.shape))
    G = as_networkx(skel, nroots=nr)
    orig_netw = nx_to_pv(G)
    splines = generate_splines(G, point_ratio=3)
    splines = [spl.interpolate(orig_netw, strategy="closest_point") for spl in splines]
    lines = [pv.lines_from_points(spl.points) for spl in splines]
    for l, spl in zip(lines, splines):
        l["radius"] = spl["radius"]
        l["root"] = spl["root"]
    smooth_netw = pv.MultiBlock(lines).combine().clean()
    orig_netw.save(f"../mesh/networks/{name}.vtk")
    smooth_netw.save(f"../mesh/networks/{name}_smooth.vtk")
    netw_tubes = as_tubes(splines).combine()
    netw_tubes.save(f"../mesh/networks/{name}_tubes.vtk")
    netw_tubes.plot(off_screen=True, screenshot=f"../plots/{name}_network.png", zoom=1.6)
    #plot_radii([s for s in splines if s.number_of_points > 4], f"../plots/{name}_arc_radii.png")
