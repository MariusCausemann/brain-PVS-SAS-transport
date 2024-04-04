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

    # Use networkx directed graph representation of mesh for convenience.
    G = nx.DiGraph()
    G.add_edges_from(mesh.cells())
    degrees = np.array([G.degree(i) for i in range(G.number_of_nodes())])
    print("Number of nodes = ", mesh.num_vertices())
    print("Number of edges = ", mesh.num_cells())

    # FIXME: Identify the network supply nodes (Currently identified
    # by a semi-visual inspection.)
    (r0, r1, r2) = [4094, 7220, 7974]

    # Compute shortest path from each of the supply nodes (r0, r1, r2)
    # to the other(s) via Dijkstra's algorithm. (Undirected graph to
    # make sure that the path is found.)
    H = G.to_undirected()
    p01 = nx.dijkstra_path(H, r0, r1)
    p02 = nx.dijkstra_path(H, r0, r2)
    p12 = nx.dijkstra_path(H, r1, r2)

    # Store paths
    print("Storing shortest paths to %s/..." % output)
    mf = df.MeshFunction("size_t", mesh, 0, 0)
    mf.array()[p01] = 1
    df.File(os.path.join(output, "path_r0r1.pvd")) << mf

    mf = df.MeshFunction("size_t", mesh, 0, 0)
    mf.array()[p02] = 2
    df.File(os.path.join(output, "path_r0r2.pvd")) << mf

    mf = df.MeshFunction("size_t", mesh, 0, 0)
    mf.array()[p12] = 3
    df.File(os.path.join(output, "path_r1r2.pvd")) << mf

if __name__ == '__main__':

    map_vasculature("../mesh/networks/arteries_smooth.vtk", "paths")
    #map_vasculature("../mesh/networks/venes_smooth.vtk")

