# mamba activate pvs_transport_env

from itertools import combinations
import os.path

import networkx as nx
import dolfin as df
import numpy as np

from solver import read_vtk_network

def compute_shortest_paths(mesh, H, roots, output="_tmp"):
    # Find and mark shortest path from each of the nodes in G to the
    # 'roots'.

    # Helper function to compute path length/weight based on edge data
    def path_length(G, path):
        length = 0.0
        for (i, _) in enumerate(path[:-1]):
            length += G.get_edge_data(path[i], path[i+1])["weight"]
        return length

    # Create a vertex function that represents whether r0, r1 or r2 is
    # the closest node for each (terminal) node
    mf = df.MeshFunction("size_t", mesh, 0, 0)

    # labels map (0, 1, ..., n_r) -> (index_0, index_1, ..., index_r)
    labels = dict(zip(range(len(roots)), roots))

    # Get all terminal nodes to begin with
    nodes = [t for t in H if H.degree(t) == 1]

    # ... and in addition all nodes in shortests paths betweeen
    # pairs of roots
    (r0, r1, r2) = roots
    pairs = combinations(roots, 2) 
    for (r0, r1) in pairs: 
        path = nx.dijkstra_path(H, r0, r1)
        nodes.extend(path)
    nodes = set(nodes)

    # Iterate over these nodes in an attempt to cover whole network
    print("Computing nearest (weighted) supply node for %d nodes " % len(nodes))
    for node in nodes:
        paths = [nx.dijkstra_path(H, node, r, weight="weight")
                 for r in roots]
        path_lengths = np.array([path_length(H, p) for p in paths])
        shortest = np.argmin(path_lengths)
        mf.array()[paths[shortest]] = labels[shortest]

    # Check that all nodes have been marked
    if not np.all(mf.array()):
        print("WARNING: Not all nodes have been marked, still zeros in mf.")

    # Store computed data as XDMF (for easy reading/writing)
    filename = os.path.join(output, "nearest_supply_nodes.xdmf")
    with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
        xdmf.write(mf)
        
    return mf


# Recursive function used to reduce full network to minimal set of
# nodes and edges
def add_branches(tree, a0, a, a_, Gr):
    degree = tree.degree(a)

    # Handle roots and leafs:
    if degree == 1:
        if a0 == a:
            (_, b) = list(tree.edges(a))[0]
            add_branches(tree, a0, b, a, Gr)
        else:
            Gr.add_edge(a0, a)
            
    # If we are on a path (degree == 2)
    if degree == 2:
        (e0, e1) = list(tree.edges(a))
        # Continue down the new edge (not including a_)
        _, b = e0 if a_ in e1 else e1
        assert (a == _), "Error in path traversal"
        add_branches(tree, a0, b, a, Gr)
        
    # Ok, at a bifurcation
    if degree == 3:
        Gr.add_edge(a0, a) 
        # Get the other edges
        (e1, e2) = [e for e in tree.edges(a) if not a_ in e]
        assert (a == e1[0] and a == e2[0]), "Error in tree traversal"
        add_branches(tree, a, e1[1], a, Gr)
        add_branches(tree, a, e2[1], a, Gr)

    # Ok, at a bifurcation with more than 4
    if degree == 4:
        Gr.add_edge(a0, a) 
        # Get the other edges
        print(tree.edges(a))
        (e1, e2, e3) = [e for e in tree.edges(a) if not a_ in e]
        assert (a == e1[0] and a == e2[0] and a == e3[0]), "Error in tree traversal"
        add_branches(tree, a, e1[1], a, Gr)
        add_branches(tree, a, e2[1], a, Gr)
        add_branches(tree, a, e3[1], a, Gr)

    if degree > 4:
        raise Exception("degree > 4")

def map_vasculature(fname, output):

    # Read network information from 'fname' -- argh, don't rescale it.
    print("Reading from %s" % fname)
    mesh, radii, roots = read_vtk_network(fname, rescale_mm2m=False)
    print("Number of network nodes = ", mesh.num_vertices())
    print("Number of network edges = ", mesh.num_cells())
    mesh.init()

    # Compute the length of each cell (lazy version)
    DG0 = df.FunctionSpace(mesh, "DG", 0)
    v = df.TestFunction(DG0)
    cell_lengths = df.assemble(1*v*df.dx(domain=mesh)).get_local()
    
    # Use networkx directed graph representation of mesh for
    # convenience. Add the size of each cell as edge weights.
    G = nx.DiGraph()
    for i, (n0, n1) in enumerate(mesh.cells()):
        G.add_edge(n0, n1, weight=cell_lengths[i])

    # Identify the network supply nodes (Currently identified by a
    # semi-visual inspection.)
    (r0, r1, r2) = [4094, 7220, 7974]

    # Easier to work with the undirected graph from here on.
    H = G.to_undirected()

    # Compute (or read pre-computed) paths to supply nodes
    #mf = compute_shortest_paths(mesh, H, (r0, r1, r2), output=output)
    filename = os.path.join(output, "nearest_supply_nodes.xdmf")
    mf = df.MeshFunction("size_t", mesh, 0, 0)
    with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
        xdmf.read(mf)

    #mf = df.MeshFunction("size_t", mesh, 0, 0)
    #mf.array()[3712] = 3712
    #mf.array()[9697] = 9697
    #file = df.File(os.path.join(output, "nodes_with_degree4.pvd"))
    #file << mf
    
    # Extract subgraphs by taking all nodes marked by r0, r1 and r2
    # separately
    for r in (r0, r1, r2):
        print("r = ", r)
        # Extract subgraph tree defined by marking
        subnodes = np.where(mf.array() == r)[0]
        tree = H.subgraph(subnodes)

        # Create reduced graph by adding bifurcations and terminals as nodes
        Gr = nx.DiGraph()
        bifurcations = [b for b in tree if tree.degree(b) >= 3]
        terminals = [t for t in tree if tree.degree(t) == 1]
        Gr.add_nodes_from(bifurcations)
        Gr.add_nodes_from(terminals)
        add_branches(tree, r, r, r, Gr)
        print("G (r = %d has %d nodes and %d edges", r, Gr.num_nodes(), Gr.num_edges())
        
if __name__ == '__main__':

    map_vasculature("../mesh/networks/arteries_smooth.vtk",
                    "../mesh/networks/common_paths")
    #map_vasculature("../mesh/networks/venes_smooth.vtk")

        
    # # Iterate over edges in subtree (depth-first)
    # a0 = r
    # for e in nx.edge_dfs(tree, source=r):
    #     (a, b) = e
    #     print("e = ", e)
    #     print("a.degree() = ", tree.degree(a))
    #     print("b.degree() = ", tree.degree(b))
    #     if tree.degree(b) == 3:
    #         b0 = b
    #         Gr.add_edge(a0, b0)
    #         a0 = b0
    #         exit()   
            
                
