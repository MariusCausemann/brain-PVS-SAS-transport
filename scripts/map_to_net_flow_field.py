import os.path

import networkx as nx
import dolfin as df
import numpy as np

from solver import read_vtk_network

def map_vasculature(fname, output):

    # Read network information from 'fname' -- argh, don't rescale it.
    print("Reading from %s" % fname)
    mesh, radii, roots = read_vtk_network(fname, rescale_mm2m=False)
    mesh.init()

    # Compute the length of each cell (lazy version)
    DG0 = df.FunctionSpace(mesh, "DG", 0)
    v = df.TestFunction(DG0)
    cell_lengths = df.assemble(1*v*df.dx(domain=mesh)).get_local()
    
    # Use networkx directed graph representation of mesh for
    # convenience. Add the size of each cell as edge weights. (MER has
    # not found an interface for doing this without adding each edge
    # separately.)
    G = nx.DiGraph()
    for i, (n0, n1) in enumerate(mesh.cells()):
        G.add_edge(n0, n1, weight=cell_lengths[i])
    
    print("Number of network nodes = ", mesh.num_vertices())
    print("Number of network edges = ", mesh.num_cells())

    # FIXME: Identify the network supply nodes (Currently identified
    # by a semi-visual inspection.)
    (r0, r1, r2) = [4094, 7220, 7974]

    # Get all terminal nodes (degree == 1)
    degree = np.array([d for (n, d) in G.degree()])
    terminals = np.where(degree == 1)[0]
    print(terminals)

    # Easier to work with the undirected graph here.
    H = G.to_undirected()
    
    # Compute shortest path from each of the terminals to the supply nodes.
    def path_length(G, path):
        length = 0.0
        for (i, v) in enumerate(path[:-1]):
            v0 = path[i]
            v1 = path[i+1]
            length += G.get_edge_data(v0, v1)["weight"]
        return length

    # Create a vertex function that holds whether r0 (1), r1 (2) or r2 (3) is the
    # closest for each terminal node 
    mf = df.MeshFunction("size_t", mesh, 0, 0)

    labels = {0: r0, 1: r1, 2: r2}

    for t in terminals:
        print("Computing nearest supply node for terminal node ", t)
        paths = [nx.dijkstra_path(H, t, r, weight="weight") for r in (r0, r1, r2)]
        path_lengths = np.array([path_length(H, p) for p in paths])
        print("path_lengths = ", path_lengths)
        shortest = np.argmin(path_lengths)
        print("shortest = ", shortest)
        mf.array()[paths[shortest]] = labels[shortest]
        
    df.File(os.path.join(output, "nearest_supply_node.pvd")) << mf

    # In addition ... 
    # Compute shortest path from each of the supply nodes (r0, r1, r2)
    # to the other(s) via Dijkstra's algorithm. (Undirected graph to
    # make sure that the path is found.)
    p01 = nx.dijkstra_path(H, r0, r1)
    p02 = nx.dijkstra_path(H, r0, r2)
    p12 = nx.dijkstra_path(H, r1, r2)

    # Find common node n0 that is in all three sets:
    n0 = set(p01).intersection(set(p02)).intersection(set(p12)).pop()
    print("Common node n0 = ", n0)
    
    p0 = nx.dijkstra_path(H, r0, n0)
    p1 = nx.dijkstra_path(H, r1, n0) 
    p2 = nx.dijkstra_path(H, r2, n0)
   
    # Store paths
    print("Storing shortest paths to %s/..." % output)
    mf = df.MeshFunction("size_t", mesh, 0, 0)
    mf.array()[p0] = r0
    mf.array()[p1] = r1
    mf.array()[p2] = r2
    df.File(os.path.join(output, "shortest_path_between_supply_nodes.pvd")) << mf

if __name__ == '__main__':

    map_vasculature("../mesh/networks/arteries_smooth.vtk", "../mesh/networks/common_paths")
    #map_vasculature("../mesh/networks/venes_smooth.vtk")

