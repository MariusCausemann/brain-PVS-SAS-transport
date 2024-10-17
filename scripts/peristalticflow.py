from itertools import combinations
import os.path
import pylab
import sys
import pyvista as pv

import networkx as nx
import dolfin as df
import numpy as np

#from solver import read_vtk_network

import pvs_network_netflow as pnf

# FIXME: Automate
SUPPLY_NODES = [8065, 7173, 4085]

def my_read_vtk_network(filename, rescale_mm2m=True):
    """Read the VTK file given by filename, return a FEniCS 1D Mesh representing the network, a FEniCS MeshFunction (double) representing the radius of each vessel segment (defined over the mesh cells), and a FEniCS MeshFunction (size_t) defining the roots of the network (defined over the mesh vertices, roots are labelled by 2 or 1.) 

    rescale_mm2m is set to True by default, in which case the information on file is rescaled by a factor 1.e-3. MER: This is error-prone, I suggest no rescaling behind the scenes.
"""
    print("Reading network mesh from %s" % filename)
    netw = pv.read(filename)
    if rescale_mm2m:
        netw.points *= 1e-3 # scale to m
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

    # Cell Function
    radii = df.MeshFunction('double', mesh, 1, 0)
    roots = df.MeshFunction('size_t', mesh, 0, 0)

    roots.array()[:] = netw["root"]
    netw = netw.point_data_to_cell_data()

    if rescale_mm2m:
        radii.array()[:] = netw["radius"] * 1e-3 # scale to m
    else:
        radii.array()[:] = netw["radius"]

    print("... with %d nodes, %d edges" % (mesh.num_vertices(), mesh.num_cells()))
    return mesh, radii, roots


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
    longest_shortest_path = 0.0
    for i in nodes:
        paths = [nx.dijkstra_path(G, i, i0, weight="length") for i0 in roots]
        path_lengths = np.array([path_length(G, p) for p in paths])
        shortest = np.argmin(path_lengths)
        longest_shortest_path = max(min(path_lengths), longest_shortest_path)
        mf.array()[paths[shortest]] = labels[shortest]
    print("Longest shortest path is %g (mm)" % longest_shortest_path)

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

def add_branches(G, a0, a, a_, T, indices, radii, lengths, index_map,
                 downstream):
    """For a given graph G with current root node a0, current node a
    and previous node a_, compute (by recursion) 
    - T: a reduced graph including 'radius' edge data , 
    - index_map: a MeshFunction map from G's cell index to T's edge index
    - downstream: a MeshFunction map over G's cells: +1 if cell is
      aligned downstream (cell tangent points from root towards
      leaves), -1 otherwise

    The given indices (list), radii (list), lengths (list) are helper
    variables to keep track of the data associated with the paths
    traversed between bifurcation points or root nodes.

    """

    # counter counts the cells in T (tracks the cell indices in T).
    global counter

    # We are now visiting node a
    # Current path is (a0, ..., a)
    # Previous edge is (a_, a)
    # When visiting node a, we take responsibility for adding previous
    # edge info to indices and radii.
    
    # Handle node a based on its graph degree
    degree = G.degree(a)

    # At the (single) root node, we just get started down the G
    if degree == 1 and a0 == a:
        (_, b) = list(G.edges(a))[0]
        add_branches(G, a0, b, a, T, [], [], [], index_map, downstream)
        return
        
    # Update indices and radii with info about the edge (a_, a) we
    # just traversed
    index = G.get_edge_data(a_, a)["index"]
    indices += [index]
    radii += [G.get_edge_data(a_, a)["radius"]]
    lengths += [G.get_edge_data(a_, a)["length"]]

    # Set downstream for this cell as 1 if aligned downstream or -1 if not
    v = index_map.mesh().cells()[index]
    if (a_ == v[0] and a == v[1]):
        downstream.array()[index] = +1.0
    elif (a_ == v[1] and a == v[0]):
        downstream.array()[index] = -1.0
    else:
        raise Exception("Something is rotten with the mesh vs G (%r)" % v)

    # Leaf node at end of branch, add new edge to minimal tree
    if degree == 1 and not (a0 == a):

        # Add new edge to minimal tree
        counter += 1
        path_radius = average_radius(radii)
        path_length = sum(lengths)
        T.add_edge(a0, a,
                   index=counter, radius=path_radius, length=path_length)
        
        # Update index map
        index_map.array()[indices] = counter
        
    # If we are on a path (degree == 2)
    if degree == 2:

        # Find which is the new edge (the edge not including a_)
        (e0, e1) = list(G.edges(a))
        _, b = e0 if a_ in e1 else e1
        assert (a == _), "Assumption error in path traversal"

        # Continue down the new edge 
        add_branches(G, a0, b, a, T, indices, radii, lengths, index_map,
                     downstream)
        
    # Ok, at a bifurcation
    if degree == 3:

        # Add new edge to minimal tree T
        counter += 1
        path_radius = average_radius(radii)
        path_length = sum(lengths)
        T.add_edge(a0, a,
                   index=counter, radius=path_radius, length=path_length)

        # Update index map
        index_map.array()[indices] = counter

        # Get the other edges 
        (e1, e2) = [e for e in G.edges(a) if not a_ in e]
        assert (a == e1[0] and a == e2[0]), "Assumption error in tree traversal"

        add_branches(G, a, e1[1], a, T, [], [], [], index_map, downstream)
        add_branches(G, a, e2[1], a, T, [], [], [], index_map, downstream)
        
    # Ok, at a bifurcation with more than 4. This shouldn't happen really
    if degree == 4:

        # Add this edge to T 
        counter += 1
        path_radius = average_radius(radii)
        path_length = sum(lengths)
        T.add_edge(a0, a,
                   index=counter, radius=path_radius, length=path_length)

        # Update index map
        index_map.array()[indices] = counter

        # Get the other edges
        (e1, e2, e3) = [e for e in G.edges(a) if not a_ in e]
        assert (a == e1[0] and a == e2[0] and a == e3[0]), \
            "Assumption error in G traversal"
        print("a ", a)
        print("e1[0] = ", e1[1])
        print("e2[0] = ", e2[1])
        print("e3[0] = ", e3[1])
        
        if (a == 8103):
            add_branches(G, a, 8201, a, T, [], [], [], index_map, downstream)
            add_branches(G, a, 12073, a, T, [], [], [], index_map, downstream)
            print("WARNING: Ignoring edge: (degree > 3) ", e2)
        else:
            add_branches(G, a, e1[1], a, T, [], [], [], index_map, downstream)
            add_branches(G, a, e2[1], a, T, [], [], [], index_map, downstream)
            print("WARNING: Ignoring edge: (degree > 3) ", e3)
        #add_branches(G, a, e3[1], a, T, [], [], [], index_map, downstream)

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

    beta = 3.0
    for G in Gs:
        (indices, paths, r_o, r_e, L) = graph_to_bifurcations(G, 0, beta)
        print("indices = ", indices)
        print("paths = ", paths)
        print("r_o = ", r_o)
        print("r_e = ", r_e)
        print("L = ", L)

def extract_minimal_tree(G, mesh, i0, downstream=None):  
    """Extract the minimal tree T from the graph G using i0 as the supply
    node, including a map M from edge indices in G/mesh to cell
    indices in T. Returns T and M.

    Use downstream map if given, otherwise create and populate new
    MeshFunction over cells.

    """

    # Start creating reduced graph by adding bifurcations and
    # terminals as nodes. FIXME: Wouldn't it be a good idea to also
    # add spatial coordinates here?
    T = nx.Graph()
    #bifurcations = [i for i in G if G.degree(i) >= 3]
    #terminals = [i for i in G if G.degree(i) == 1]
    #T.add_nodes_from(bifurcations)
    #T.add_nodes_from(terminals)

    # ... and then adding edges to the reduced graph
    # ... making sure to also make a map from cell indices in the
    # original tree to an index in the new graph

    UNDEFINED = mesh.num_cells()
    cell_index_map = df.MeshFunction("size_t", mesh, 1, UNDEFINED)
    if not downstream:
        downstream = df.MeshFunction("double", mesh, 1, 0)
    global counter
    counter = -1

    add_branches(G, i0, i0, i0, T, [], [], [], cell_index_map, downstream)
    nv = T.number_of_nodes()
    ne = T.number_of_edges()
    print("... extracted minimal tree T with %d nodes and %d edges" %
          (T.number_of_nodes(), T.number_of_edges()))

    assert (nv == (ne + 1)), "Number of nodes and edges do not match!"
    
    return T, cell_index_map, downstream

def compute_subtrees(filename, output):

    # Label the network supply nodes (FIXME: automate)
    supply_nodes = SUPPLY_NODES
    
    # Read network information from file. Never rescale stuff behind
    # the scenes.
    mesh, radii, _ = my_read_vtk_network(filename, rescale_mm2m=False)
    mesh.init()

    # Convert mesh to weighted graph
    G = mesh_to_weighted_graph(mesh, radii)

    # Compute and store shortest paths. If read = True, just read from
    # file (assuimng that it is precomputed)
    mf = mark_shortest_paths(mesh, G, supply_nodes, output)

    # Extract subtrees from G corresponding for each supply node, and
    # compute minimal tree representation for each
    downstream = df.MeshFunction("double", mesh, 1, 0)
    for i0 in supply_nodes:

        print("Computing minimal subtree starting at %d" % i0)
        subnodes = np.where(mf.array() == i0)[0]
        G0 = G.subgraph(subnodes).copy()
        T, cell_index_map, _ = extract_minimal_tree(G0, mesh, i0, downstream)
        
        # Store the new graph 
        filename = os.path.join(output, "minimal_tree_%d.graphml" % i0)
        print("... storing minimal tree T representation to %s" % filename)
        nx.write_graphml(T, filename, infer_numeric_types=True)
                                
        # ... and the G-to-T cell index map
        filename = os.path.join(output, "original_to_minimal_map_%d.xdmf" % i0)
        with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
            print("... storing index map (mesh -> T) to %s" % filename)
            xdmf.write(cell_index_map)

    # ... and the downstream orientation map: cell orientation vs
    # downstream orientation (aligned -> 1, not aligned = -1)
    filename = os.path.join(output, "downstream_map.xdmf")
    with df.XDMFFile(mesh.mpi_comm(), filename) as xdmf:
        print("... storing downstream map (mesh -> +-1) to %s" % filename)
        xdmf.write(downstream)
        
def dg0_to_mf(u):
    mesh = u.function_space().mesh()
    mf = df.MeshFunction("double", mesh, mesh.topology().dim(), 0)
    mf.array()[:] = u.vector().get_local()
    return mf

def print_stats(name, L, unit):
    print("%s (avg, std, min, max): %.4g, %.4g, %.4g, %.4g (%s)" %
          (name, np.average(L), np.std(L), L.min(), L.max(), unit))

def compute_pvs_flow(meshfile, output, args):

    # Label the network supply nodes (FIXME: automate)
    roots = SUPPLY_NODES

    # Read mesh from file. Never rescale stuff behind the scenes.
    mesh, radii, _ = my_read_vtk_network(meshfile, rescale_mm2m=False)
    mesh.init()
    DG0 = df.FunctionSpace(mesh, "DG", 0)

    print_stats("Vascular radii r_o", radii.array(), "mm")
    
    # Note that the asymptotic estimate is is derived under the
    # assumption that when k L = 2 pi/lmbda L = O(1) i.e. when lmbda ~
    # 2 pi L.

    # Specify the relative PVS width and other parameters
    beta = args.beta             # outer radius = beta*(inner radius)
    f = args.frequency           # Frequency
    lmbda = args.wavelength      # Wave length (mm)
    varepsilon = args.amplitude  # Relative wave amplitude

    omega = 2*np.pi*f      # Angular frequency (Hz)
    k = 2*np.pi/lmbda      # Wave number (1/mm)
    
    Q = df.Function(DG0)
    u = df.Function(DG0)
    Ls = []

    # Short-hand to compute annular cross-section area for cell i
    area = lambda i: np.pi*(beta**2-1.0)*radii[i]**2 

    print("Computing time-average perivascular flow rates...")
    for i0 in roots:
        # Read minimal subtree from file. Note that graphml converts
        # ints to strings ...
        graphfile = os.path.join(output, "minimal_tree_%d.graphml" % i0)
        T = nx.read_graphml(graphfile)
        print("\tT (%d)" % i0, T)
        
        # Map minimal subtree T into PVS net flow data representation
        (indices, paths, r_o, r_e, L) = graph_to_bifurcations(T, i0, beta)
        Ls += L
        
        # Compute the PVS net flow rates in T
        network_data = (indices, paths, r_o, r_e, L, k, omega, varepsilon)
        avg_Q, avg_u = pnf.estimate_net_flow(network_data)
        
        # Read cell_index_map from file for mapping back to the original mesh 
        mapfile = os.path.join(output, "original_to_minimal_map_%d.xdmf" % i0)
        index_map = df.MeshFunction("size_t", mesh, 1, 0)
        with df.XDMFFile(mesh.mpi_comm(), mapfile) as xdmf:
            print("\t... reading index map (mesh -> T) from %s" % mapfile)
            xdmf.read(index_map)

        # Transfer the computed flow rates back to the mesh. Compute
        # velocities by dividing by flux by cross-section area
        # cell-wise.
        UNDEFINED = mesh.num_cells()
        for i, _ in enumerate(mesh.cells()):
            T_i = index_map[i]
            if (T_i < UNDEFINED):
                avg_Q_i = avg_Q[index_map[i]]
                Q.vector()[i] = avg_Q_i
                u.vector()[i] = avg_Q_i/area(i)
                
    print_stats("Vascular branch lengths L", np.array(Ls), "mm")            
    print_stats("... k L", k*np.array(Ls), "AU, (k = %3.4g)" % k )
    print_stats("<Q'_i>", Q.vector(), "mm^3/s")
    print_stats("<u'_i>", u.vector(), "mm/s")
    print_stats("<u'_i>", u.vector()*1.e3, "mum/s")

    fluxfile = os.path.join(output, "pvs_Q.xdmf")
    with df.XDMFFile(mesh.mpi_comm(), fluxfile) as xdmf:
        print("Saving net flux to %s" % fluxfile)
        xdmf.write_checkpoint(Q, "Q", 0.0)

    fluxfile = os.path.join(output, "pvs_Q_mf.xdmf")
    with df.XDMFFile(mesh.mpi_comm(), fluxfile) as xdmf:
        print("Saving net flux (as mf) to %s" % fluxfile)
        xdmf.write(dg0_to_mf(Q))
        
    ufile = os.path.join(output, "pvs_u.xdmf")
    with df.XDMFFile(mesh.mpi_comm(), ufile) as xdmf:
        print("Saving net velocity to %s" % ufile)
        xdmf.write_checkpoint(u, "u", 0.0)

    ufile = os.path.join(output, "pvs_u_mf.xdmf")
    with df.XDMFFile(mesh.mpi_comm(), ufile) as xdmf:
        print("Saving net velocity (as mf) to %s" % ufile)
        xdmf.write(dg0_to_mf(u))
                
def run_all_tests():
    test_graph_to_bifurcations()

def main(args):

    # Define input network (.vtk) and output directory
    filename = "../mesh/networks/arteries_smooth.vtk"
    output = "../mesh/networks/arterial_trees"

    # Only need to compute subtree information once for each mesh
    if args.recompute:
        compute_subtrees(filename, output)
        print("")
        
    #compute_pvs_flow(filename, output, args)
    
if __name__ == '__main__':

    # TODO: Make and automatically run tests that checks these
    # functions on a simple graphs/meshes/networks and on our favorite
    # image-based network

    # To run:
    # $ mambda .... (FIXME) marie_environment.xml
    # $ mamba activate pvs_transport_env
    # $ pip install git+https://github.com/MiroK/fenics_ii.git@ufl_legacy
    # 
    # $ python3 peristalticflow.py --frequency 1.0 --wavelength 2000.0 --amplitude 0.01 --beta 3.0 --recompute # Cardiac
    # $ python3 peristalticflow.py --frequency 0.1 --wavelength 80.0 --amplitude 0.1 --beta 3.0 
    #
    # --recompute is only needed upon first run)
    
    import argparse

    parser = argparse.ArgumentParser(description="Compute estimate of net flow in a perivascular network", epilog="Run with --recompute on first go.")
    parser.add_argument('--recompute', action="store_true")
    parser.add_argument('--run-tests', action="store_true", help="Just run basic test.")
    parser.add_argument('--frequency', action="store", default=1.0, help="Vascular wave frequency (Hz)", type=float)
    parser.add_argument('--wavelength', action="store", default=2000.0, help="Vascular wave wavelength (mm)", type=float)
    parser.add_argument('--amplitude', action="store", default=0.01, help="Vascular wave relative amplitude (relative to (inner) vascular radius)", type=float)
    parser.add_argument('--beta', action="store", default=3.0, help="Ratio outer-to-inner vessel radio (PVS width + 1)", type=float)
    parser.add_argument('--tag', action="store", default="tmp", help="Tag")

    args = parser.parse_args()
    if args.run_tests:
        run_all_tests()
        exit()
    main(args)

