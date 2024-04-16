# mamba activate pvs_transport_env

from itertools import combinations
import os.path

import networkx as nx
import dolfin as df
import numpy as np

from solver import read_vtk_network

import pvs_network_netflow

def mark_shortest_paths(mesh, G, roots, output, use_precomputed=False):
    """Identify the shortest path (weighted by "length") from each node in
    G to a root node (in roots). Create a node-based mesh function mf
    that labels each note by its nearest root. Returns mf.
    """

    # Helper function to compute path length based on edge data
    def path_length(G, path):
        length = 0.0
        for (i, _) in enumerate(path[:-1]):
            length += G.get_edge_data(path[i], path[i+1])["length"]
        return length

    mf = df.MeshFunction("size_t", mesh, 0, 0)
    filename = os.path.join(output, "nearest_supply_nodes.xdmf")

    # Read from file if relevant
    if use_precomputed:
        print("Reading precomputed nearest root labels from %s" % filename)
        with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
            xdmf.read(mf)
        return mf

    # labels map (0, 1, ..., n_r) -> (index_0, index_1, ..., index_r)
    labels = dict(zip(range(len(roots)), roots))

    # Get all terminal nodes to begin with
    nodes = [t for t in G if G.degree(t) == 1]

    # ... and in addition all nodes in shortests paths betweeen
    # pairs of roots, to cover whole network
    pairs = combinations(roots, 2) 
    for (i0, i1) in pairs: 
        path = nx.dijkstra_path(G, i0, i1)
        nodes.extend(path)
    nodes = set(nodes)

    # For these nodes, identify which root node that is the closest by
    # computing the shortest shortest path, and label all nodes in
    # path by that root
    print("Finding nearest root for %d nodes, this may take some time..." % len(nodes))
    for i in nodes:
        paths = [nx.dijkstra_path(G, i, i0, weight="length") for i0 in roots]
        path_lengths = np.array([path_length(G, p) for p in paths])
        shortest = np.argmin(path_lengths)
        mf.array()[paths[shortest]] = labels[shortest]

    # Check that all nodes have been marked
    if not np.all(mf.array()):
        print("WARNING: Not all nodes have been marked, still zeros in mf.")

    # Store computed data as XDMF (for easy reading/writing)
    print("... saving nearest root labels as MeshFunction to %s" % filename)
    with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
        xdmf.write(mf)
        
    return mf

def average_radius(radii):
    "How to compute radius of collapsed branch? Here, we take the average."
    return np.average(np.array(radii))

def add_branches(tree, a0, a, a_, Gr, indices, radii, lengths, index_map):
    """For a given graph tree with current root node a0, current node a
    and previous node a_, compute (by recursion) 
    - Gr: a reduced graph including 'radius' edge data , 
    - index_map: a MeshFunction map from tree's cell index to Gr's cell index

    The given indices (list), radii (list), lengths (list) are helper
    variables to keep track of the data associated with the paths
    traversed between bifurcation points or root nodes.
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
        add_branches(tree, a0, b, a, Gr, [], [], [], index_map)
        return
        
    # Update indices and radii with info about the edge (a_, a) we
    # just traversed
    indices += [tree.get_edge_data(a_, a)["index"]]
    radii += [tree.get_edge_data(a_, a)["radius"]]
    lengths += [tree.get_edge_data(a_, a)["length"]]
    
    # Leaf node at end of branch, add new edge to minimal tree
    if degree == 1 and not (a0 == a):

        # Add new edge to minimal tree
        counter += 1
        path_radius = average_radius(radii)
        path_length = sum(lengths)
        Gr.add_edge(a0, a,
                    index=counter, radius=path_radius, length=path_length)
        
        # Update index map
        index_map.array()[indices] = counter

    # If we are on a path (degree == 2)
    if degree == 2:

        # Find which is the new edge (the edge not including a_)
        (e0, e1) = list(tree.edges(a))
        _, b = e0 if a_ in e1 else e1
        assert (a == _), "Assumption error in path traversal"

        # Continue down the new edge 
        add_branches(tree, a0, b, a, Gr, indices, radii, lengths, index_map)
        
    # Ok, at a bifurcation
    if degree == 3:

        # Add new edge to minimal tree Gr
        counter += 1
        path_radius = average_radius(radii)
        path_length = sum(lengths)
        Gr.add_edge(a0, a,
                    index=counter, radius=path_radius, length=path_length)

        # Update index map
        index_map.array()[indices] = counter
            
        # Get the other edges 
        (e1, e2) = [e for e in tree.edges(a) if not a_ in e]
        assert (a == e1[0] and a == e2[0]), "Assumption error in tree traversal"

        add_branches(tree, a, e1[1], a, Gr, [], [], [], index_map)
        add_branches(tree, a, e2[1], a, Gr, [], [], [], index_map)

    # Ok, at a bifurcation with more than 4
    if degree == 4:

        # Add this edge to Gr, 
        counter += 1
        path_radius = average_radius(radii)
        path_length = sum(lengths)
        Gr.add_edge(a0, a,
                    index=counter, radius=path_radius, length=path_length)

        # Update index map
        index_map.array()[indices] = counter

        # Get the other edges
        print("Degree 4 edges: ", tree.edges(a))
        (e1, e2, e3) = [e for e in tree.edges(a) if not a_ in e]
        assert (a == e1[0] and a == e2[0] and a == e3[0]), \
            "Assumption error in tree traversal"

        add_branches(tree, a, e1[1], a, Gr, [], [], [], index_map)
        add_branches(tree, a, e2[1], a, Gr, [], [], [], index_map)
        add_branches(tree, a, e3[1], a, Gr, [], [], [], index_map)

    if degree > 4:
        raise Exception("degree > 4")

def mesh_to_weighted_graph(mesh, radii):
    """Given a FEniCS Mesh mesh and an indexable array radii, create a networkx
graph G, with mesh cell sizes, indices and the corresponding
radii as the edge "length", "index" and "radius". Returns G.""" 

    # Compute the length of each cell (lazy version)
    DG0 = df.FunctionSpace(mesh, "DG", 0)
    v = df.TestFunction(DG0)
    sizes = df.assemble(1*v*df.dx(domain=mesh)).get_local()

    # Create graph and add attributes to edges
    G = nx.Graph()
    for i, (n0, n1) in enumerate(mesh.cells()):
        G.add_edge(n0, n1, length=sizes[i], index=i, radius=radii[i])

    return G
        
def graph_to_bifurcations(G, i0, relative_pvs_width):
    """Map a bifurcating tree graph G (networkx) with root node i0 into
    the list-based data structure used to compute net PVS flow by
    peristalsis.
    """

    i0 = str(i0) # FIXME: Hack due to how graphml write/reads nodes.
    
    # The bifurcation points are the degree 3 nodes
    bifurcations = [i for i in G if G.degree(i) == 3]

    # Define tuples (e0, e1, e2) where e0 and e1, e2 are the respecive
    # indices of the mother and two daughter edges at the bifurcation.
    indices = [tuple(G[i][j]["index"] for j in G.neighbors(i))
               for i in bifurcations]

    print(indices)
    exit()
    
    # Extract paths (lists of edge indices) from the root node to each
    # leaf node
    leaves = [i for i in G if G.degree(i) == 1]
    leaves.remove(i0)
    paths = []
    for l in leaves:
        crumbs = nx.dijkstra_path(G, i0, l)
        paths.append(tuple(G[i][j]["index"]
                           for (i, j) in zip(crumbs[:-1], crumbs[1:])))

    # List inner (r_o, r_1) and outer (r_e, r_2) PVS radii by edge index
    n = G.number_of_edges()
    r_o = [0,]*n 
    r_e = [0,]*n
    L = [0,]*n
    for (u, v, d) in G.edges(data=True):
        e0 = d["index"]
        r_o[e0] = d["radius"]
        r_e[e0] = relative_pvs_width*d["radius"]
        L[e0] = d["length"]
    
    return (indices, paths, r_o, r_e, L)

def test_graph_to_bifurcations():
    # Simple bifurcation
    G = nx.Graph()
    G.add_edges_from([("0", "1", {"radius": 0.1, "length": 2.0, "index": 0}),
                      ("1", "2", {"radius": 0.1, "length": 1.0, "index": 1}),
                      ("1", "3", {"radius": 0.1, "length": 1.0, "index": 2})])
    Gs = [G]

    # Simple asymmetric bifurcating tree
    G = nx.Graph()
    G.add_edges_from([("0", "1", {"radius": 0.1, "length": 2.0, "index": 0}),
                      ("1", "2", {"radius": 0.1, "length": 1.0, "index": 1}),
                      ("1", "3", {"radius": 0.1, "length": 1.0, "index": 2}),
                      ("3", "4", {"radius": 0.05, "length": 0.5, "index": 3}),
                      ("3", "5", {"radius": 0.05, "length": 0.5, "index": 4})])
    Gs += [G]

    beta = 2
    for G in Gs:
        (indices, paths, r_o, r_e, L) = graph_to_bifurcations(G, 0, beta)
        print("indices = ", indices)
        print("paths = ", paths)
        print("r_o = ", r_o)
        print("r_e = ", r_e)
        print("L = ", L)

def extract_minimal_tree(G, mesh, i0):  
    """Extract the minimal tree T from the graph G using i0 as the supply
    node, including a map M from edge indices in G/mesh to cell
    indices in T. Returns T and M.
    """

    # Start creating reduced graph by adding bifurcations and
    # terminals as nodes. FIXME: Wouldn't it be a good idea to also
    # add spatial coordinates here?
    T = nx.Graph()
    bifurcations = [i for i in G if G.degree(i) >= 3]
    terminals = [i for i in G if G.degree(i) == 1]
    T.add_nodes_from(bifurcations)
    T.add_nodes_from(terminals)

    # ... and then adding edges to the reduced graph
    # ... making sure to also make a map from cell indices in the
    # original tree to an index in the new graph
    cell_index_map = df.MeshFunction("size_t", mesh, 1, 0)
    global counter
    counter = 0
    add_branches(G, i0, i0, i0, T, [], [], [], cell_index_map)
    nv = T.number_of_nodes()
    ne = T.number_of_edges()
    print("... extracted minimal tree T with %d nodes and %d edges" % (nv, ne))

    return T, cell_index_map

def compute_subtrees(filename, output):

    # Label the network supply nodes (FIXME: automate)
    supply_nodes = [4094, 7220, 7974]

    # Read network information from file. Never rescale stuff behind
    # the scenes.
    mesh, radii, _ = read_vtk_network(filename, rescale_mm2m=False)
    mesh.init()

    # Convert mesh to weighted graph
    G = mesh_to_weighted_graph(mesh, radii)

    # Compute and store shortest paths. If read = True, just read from
    # file (assuimng that it is precomputed)
    mf = mark_shortest_paths(mesh, G, supply_nodes, output)

    # Extract subtrees from G corresponding for each supply node, and
    # compute minimal tree representation for each
    for i0 in supply_nodes:

        print("Computing minimal subtree starting at %d" % i0)
        subnodes = np.where(mf.array() == i0)[0]
        G0 = G.subgraph(subnodes).copy()
        T, cell_index_map = extract_minimal_tree(G0, mesh, i0)
        
        # Store the new graph and cell_index_map
        filename = os.path.join(output, "minimal_tree_%d.graphml" % i0)
        print("... storing minimal tree T representation to %s" % filename)
        nx.write_graphml(T, filename, infer_numeric_types=True)
                                
        # Store the cell index to reduced cell index map
        filename = os.path.join(output, "original_to_minimal_map_%d.xdmf" % i0)
        with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
            print("... storing index map (mesh -> T) to %s" % filename)
            xdmf.write(cell_index_map)

def compute_pvs_flow(meshfile, output):

    # Label the network supply nodes (FIXME: automate)
    roots = [4094, 7220, 7974]

    # Read mesh from file. Never rescale stuff behind the scenes.
    mesh, radii, _ = read_vtk_network(meshfile, rescale_mm2m=False)
    mesh.init()

    # Specify the relative PVS width
    beta = 2
    
    for i0 in roots:
        # Read minimal subtree from file. Note that graphml converts
        # ints to strings ...
        graphfile = os.path.join(output, "minimal_tree_%d.graphml" % i0)
        T = nx.read_graphml(graphfile)
            
        # Map reduced graph Gr into PVS net flow data representation
        (indices, paths, r_o, r_e, L) = graph_to_bifurcations(T, i0, beta)

def run_all_tests():
    test_graph_to_bifurcations()

def main():

    # Define input network (.vtk) and output directory
    filename = "../mesh/networks/arteries_smooth.vtk"
    output = "../mesh/networks/arterial_trees"

    # Only need to compute subtree information once for each mesh
    if False:
        compute_subtrees(filename, output)

    compute_pvs_flow(filename, output)
    
if __name__ == '__main__':

    # TODO: Make and automatically run tests that checks these
    # functions on a simple graphs/meshes/networks and on our favorite
    # image-based network

    if False:
        run_all_tests()

    if True:
        main()

