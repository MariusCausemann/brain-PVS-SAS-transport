import nibabel
import pyvista as pv
import numpy as np
import kimimaro
import os
import networkx as nx
from cloudvolume import Bbox
import pytetwild

mm2m = 1e-3

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
                             teasar_params={"scale": scale, "const": const})
    joined_skels = kimimaro.join_close_components(skels.values(), radius=50)
    ds_skel = joined_skels.downsample(2)
    return ds_skel 

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

def remove_duplicate_cells(netw):
    cells = np.array(netw.cells.reshape(-1, 3))[:,1:]
    cells.sort(axis=1)
    unique_cells = np.unique(cells, axis=0)
    netw.cells = np.pad(unique_cells, pad_width=((0,0), (1,0)), constant_values=2)

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
        t = spl.tube(scalars="radius", absolute=True)
        tubes.append(t)
    return tubes

def pvnetwork_to_polydata(netw):
    poly = pv.PolyData(netw.points, lines=netw.cells)
    for k, data in netw.point_data.items():
        poly[k] = data
    return poly

def pv_to_df(netw):
    import dolfin as df
    mesh = df.Mesh()
    ed = df.MeshEditor()
    ed.open(mesh, 'interval', 1, 3)
    ed.init_vertices(netw.number_of_points)
    cells = netw.cells.reshape((-1,3))
    ed.init_cells(cells.shape[0])

    for vid, v in enumerate(netw.points):
        ed.add_vertex(vid, v)
    for cid, c in enumerate(cells[:,1:]):
        ed.add_cell(cid, c)
    ed.close()
    return mesh

def seg_to_spline(seg, point_ratio):
    degrees = np.array([len(seg.point_neighbors(i)) for i in range(seg.n_points)])
    roots = np.where(degrees==1)[0]
    previous = None
    current = roots[0]
    pointlist = [current]
    while True:
        nei = seg.point_neighbors(current)
        if previous is not None:
            nei.remove(previous)
        if len(nei)==0:
            break
        previous = current
        current = nei[0]
        pointlist.append(current)
    spl = pv.Spline(points=seg.points[pointlist], n_points=point_ratio*seg.n_points)
    return spl.interpolate(seg, strategy="closest_point", radius=4e-3)

def get_splines(netw, point_ratio=1):
    from branch_marking import color_branches
    segments, segids, _ = color_branches(pv_to_df(netw))
    pvsegments = [netw.extract_cells(segments.array()[:]==si) for si in segids]
    return [seg_to_spline(pvs, point_ratio) for pvs in pvsegments]

def get_tubes(netw, solidify=False, **kwargs):
    tubes = as_tubes(get_splines(netw)).combine()
    if not solidify: return tubes
    return pytetwild.tetrahedralize_pv(tubes.extract_surface(), **kwargs)


def smooth_radius(netw):
    r = netw["radius"]
    rnew = np.ones_like(r)
    for i in range(netw.n_points):
        ns = netw.point_neighbors(i)
        if len(ns)<=2:
            rnew[i] = 0.5*r[i] + 0.5*(r[ns].mean())
        else:
            rnew[i] = max(r[ns + [i]])
    netw["radius"] = rnew

def smooth_coords(netw):
    coords = netw.points
    coordsnew = np.copy(coords)
    for i in range(netw.n_points):
        ns = netw.point_neighbors(i)
        # only smooth inner points, leave end and bifurcation points
        if len(ns)==2: 
            coordsnew[i,:] = 0.5*coords[i,:] + 0.5*(coords[ns,:].mean(axis=0))
    netw.points[:,:] = coordsnew

scale = 1.5
const = 3.0


if __name__=="__main__":
    os.makedirs("mesh/networks", exist_ok=True)
    os.makedirs("plots", exist_ok=True)

    files = ["data/pcbi.1007073.s007.nii.gz", "data/pcbi.1007073.s008.nii.gz"]
    names = ["arteries", "venes"]
    nroots = [3, 50]

    for file, name, nr in zip(files, names, nroots):
        data = nibabel.load(file)
        resolution = data.header["pixdim"][1:4]
        print(f"resolution: {resolution}")
        img = data.get_fdata().astype(int)
        skel = skeletonize(img, resolution)
        if name=="arteries":
            skel = skel.crop(Bbox([0, 0, 37], img.shape))
        G = as_networkx(skel, nroots=nr)
        orig_netw = nx_to_pv(G)
        orig_netw["radius"] *= mm2m
        orig_netw.scale(mm2m, inplace=True)
        remove_duplicate_cells(orig_netw)
        orig_netw.save(f"mesh/networks/{name}.vtk")
        splines = get_splines(orig_netw, point_ratio=3)
        lines = [pv.lines_from_points(spl.points) for spl in splines]
        for l, spl in zip(lines, splines):
            l["radius"] = spl["radius"]
            l["root"] = spl["root"]
        mean_radius = [l["radius"].mean() for l in lines]
        order_by_radius = np.flip(np.argsort(mean_radius))
        smooth_netw = pv.merge([lines[i] for i in order_by_radius],
                                merge_points=True).cast_to_unstructured_grid()
        for i in range(3):
            smooth_radius(smooth_netw)
        for i in range(3):
            smooth_coords(smooth_netw)
        smooth_netw.save(f"mesh/networks/{name}_smooth.vtk")
        netw_tubes = get_tubes(smooth_netw)
        netw_tubes.save(f"mesh/networks/{name}_tubes.vtk")
        print(smooth_netw.n_points)
        #netw_tubes.plot(off_screen=True, screenshot=f"plots/{name}_network.png", zoom=1.6)
        #plot_radii([s for s in splines if s.number_of_points > 4], f"../plots/{name}_arc_radii.png")

