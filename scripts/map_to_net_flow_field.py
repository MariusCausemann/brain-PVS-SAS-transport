# mamba activate pvs_transport_env

from itertools import combinations
import os.path

import networkx as nx
import dolfin as df
import numpy as np

from solver import read_vtk_network

def compute_shortest_paths(mesh, G, roots, output="_tmp"):
    # Find and mark shortest path from each of the nodes in G to the
    # 'roots'.

    # Helper function to compute path length based on edge data
    def path_length(G, path):
        length = 0.0
        for (i, _) in enumerate(path[:-1]):
            length += G.get_edge_data(path[i], path[i+1])["length"]
        return length

    # Create a vertex function that represents whether r0, r1 or r2 is
    # the closest node for each (terminal) node
    mf = df.MeshFunction("size_t", mesh, 0, 0)

    # labels map (0, 1, ..., n_r) -> (index_0, index_1, ..., index_r)
    labels = dict(zip(range(len(roots)), roots))

    # Get all terminal nodes to begin with
    nodes = [t for t in G if G.degree(t) == 1]

    # ... and in addition all nodes in shortests paths betweeen
    # pairs of roots
    pairs = combinations(roots, 2) 
    for (r0, r1) in pairs: 
        path = nx.dijkstra_path(G, r0, r1)
        nodes.extend(path)
    nodes = set(nodes)

    # Iterate over these nodes in an attempt to cover whole network
    print("Computing nearest supply node for %d nodes, this may take some time..." % len(nodes))
    for node in nodes:
        paths = [nx.dijkstra_path(G, node, r, weight="length")
                 for r in roots]
        path_lengths = np.array([path_length(G, p) for p in paths])
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

def average_radius(radii):
    "How to compute radius of collapsed branch? Here, we take the average."
    return np.average(np.array(radii))

# Recursive function used to reduce full network to minimal set of
# nodes and edges
def add_branches(tree, a0, a, a_, Gr, indices, radii, index_map):
    """For a given graph tree with current root node a0, current node a
    and previous node a_, compute (by recursion) 
    - Gr: a reduced graph including 'radius' edge data , 
    - index_map: a MeshFunction map from tree's cell index to Gr's cell index

    The given indices (list), radii (list) are helper variables to
    keep track of the data associated with the paths traversed between
    bifurcation points or root nodes.
    """

    # counter counts the cells in Gr (tracks the cell indices in Gr).
    global counter

    # We are now visiting node a
    # Current path is (a0, ..., a)
    # Previous edge is (a_, a)
    # When visiting node a, we take responsibility for adding previous
    # edge info to indices and radii.
    
    # Handle node a based on its graph degree
    degree = tree.degree(a)

    # At the (single) root node, we just get started down the tree
    if degree == 1 and a0 == a:
        (_, b) = list(tree.edges(a))[0]
        add_branches(tree, a0, b, a, Gr, [], [], index_map)
        return
        
    # Update indices and radii with info about the edge (a_, a) we
    # just traversed
    indices += [tree.get_edge_data(a_, a)["index"]]
    radii += [tree.get_edge_data(a_, a)["radius"]]
    
    # Leaf node at end of branch, add new edge to minimal tree
    if degree == 1 and not (a0 == a):

        # Add new edge to minimal tree
        counter += 1
        path_radius = average_radius(radii)
        Gr.add_edge(a0, a, index=counter, radius=path_radius)

        # Update index map
        index_map.array()[indices] = counter

    # If we are on a path (degree == 2)
    if degree == 2:

        # Find which is the new edge (the edge not including a_)
        (e0, e1) = list(tree.edges(a))
        _, b = e0 if a_ in e1 else e1
        assert (a == _), "Assumption error in path traversal"

        # Continue down the new edge 
        add_branches(tree, a0, b, a, Gr, indices, radii, index_map)
        
    # Ok, at a bifurcation
    if degree == 3:

        # Add new edge to minimal tree Gr
        counter += 1
        path_radius = average_radius(radii)
        Gr.add_edge(a0, a, index=counter, radius=path_radius)

        # Update index map
        index_map.array()[indices] = counter
            
        # Get the other edges 
        (e1, e2) = [e for e in tree.edges(a) if not a_ in e]
        assert (a == e1[0] and a == e2[0]), "Assumption error in tree traversal"

        add_branches(tree, a, e1[1], a, Gr, [], [], index_map)
        add_branches(tree, a, e2[1], a, Gr, [], [], index_map)

    # Ok, at a bifurcation with more than 4
    if degree == 4:

        # Add this edge to Gr, 
        counter += 1
        path_radius = average_radius(radii)
        Gr.add_edge(a0, a, index=counter, radius=path_radius)

        # Update index map
        index_map.array()[indices] = counter

        # Get the other edges
        print("Degree 4 edges: ", tree.edges(a))
        (e1, e2, e3) = [e for e in tree.edges(a) if not a_ in e]
        assert (a == e1[0] and a == e2[0] and a == e3[0]), \
            "Assumption error in tree traversal"

        add_branches(tree, a, e1[1], a, Gr, [], [], index_map)
        add_branches(tree, a, e2[1], a, Gr, [], [], index_map)
        add_branches(tree, a, e3[1], a, Gr, [], [], index_map)

    if degree > 4:
        raise Exception("degree > 4")

def map_vasculature(fname, output, supply_nodes, compute_paths=True):

    # Read network information from 'fname' -- argh, don't rescale
    # it. Never rescale stuff behind the scenes.
    print("Reading from %s" % fname)
    mesh, radii, roots = read_vtk_network(fname, rescale_mm2m=False)
    print("Number of network nodes = ", mesh.num_vertices())
    print("Number of network edges = ", mesh.num_cells())
    mesh.init()

    # Compute the length of each cell (lazy version)
    DG0 = df.FunctionSpace(mesh, "DG", 0)
    v = df.TestFunction(DG0)
    cell_lengths = df.assemble(1*v*df.dx(domain=mesh)).get_local()
    
    # Use networkx graph representation of mesh for convenience. Add
    # the size of each cell as edge length, and add cell index and
    # inner radius as auxiliary edge data.
    G = nx.Graph()
    for i, (n0, n1) in enumerate(mesh.cells()):
        G.add_edge(n0, n1, length=cell_lengths[i], index=i, radius=radii[i])
    
    # Compute or read pre-computed paths to supply nodes
    filename = os.path.join(output, "nearest_supply_nodes.xdmf")
    if compute_paths:
        mf = compute_shortest_paths(mesh, G, supply_nodes, output=output)
        print("Storing to %s" % filename)
        with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
            xdmf.write(mf)
    else:
        print("Reading from %s" % filename)
        mf = df.MeshFunction("size_t", mesh, 0, 0)
        with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
            xdmf.read(mf)

    # Extract subgraphs by taking all nodes marked by each supply node
    # separately
    for r in supply_nodes:
        print("Extracting reduced graphs with supply node r = ", r)

        # Extract subgraph tree defined by current marker r
        subnodes = np.where(mf.array() == r)[0]
        tree = G.subgraph(subnodes).copy()
        
        # Create directed reduced graph by adding bifurcations and
        # terminals as nodes. FIXME: Wouldn't it be a good idea to also add
        # spatial coordinates here?
        Gr = nx.DiGraph()
        bifurcations = [b for b in tree if tree.degree(b) >= 3]
        terminals = [t for t in tree if tree.degree(t) == 1]
        Gr.add_nodes_from(bifurcations)
        Gr.add_nodes_from(terminals)

        # ... and then adding edges to the reduced graph
        # ... making sure to also make a map from cell indices in the
        # original tree to an index in the new graph
        cell_index_map= df.MeshFunction("size_t", mesh, 1, 0)
        global counter
        counter = 0
        add_branches(tree, r, r, r, Gr, [], [], cell_index_map)
        print("G_r (r = %d) has %d nodes and %d edges" % 
              (r, Gr.number_of_nodes(), Gr.number_of_edges()))

        # Store the new graph
        filename = os.path.join(output, "reduced_graph_%d.graphml" % r)
        print("Storing reduced graph Gr to %s" % filename)
        nx.write_graphml(Gr, filename)
                                
        # Store the cell index to reduced cell index map
        filename = os.path.join(output, "reduced_index_map_%d.xdmf" % r)
        with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
            print("Storing index map (tree -> Gr) to %s" % filename)
            xdmf.write(cell_index_map)

def map_tree_into_bifurcations():
    pass
            
def main(compute_paths=False, compute_tree=False):

    # Compute reduced representation of arterial tree
    network = "../mesh/networks/arteries_smooth.vtk"
    output = "../mesh/networks/arterial_trees"

    # Identify the network supply nodes (Currently identified by a
    # semi-visual inspection.)
    supply_nodes = [4094, 7220, 7974]

    # Map arterial network into a set of bifurcating trees, and store
    if compute_tree:
        map_vasculature(network, output, supply_nodes, compute_paths)
    
    for r in supply_nodes:
        filename = os.path.join(output, "reduced_graph_%d.graphml" % r)
        Gr = nx.read_graphml(filename)
        print("Reading in", Gr, "from", filename)

    # Map reduced graph Gr into PVS net flow data representation

    
        
    
if __name__ == '__main__':

    # TODO: Make a few tests that checks these functions on a simple
    # graphs/meshes/networks
    
    # Set these to True when running the first time.
    main(compute_paths=False, compute_tree=False)
    #map_vasculature("../mesh/networks/venes_smooth.vtk")

        
