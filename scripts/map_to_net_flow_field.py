import networkx as nx
import dolfin
import numpy as np

from solver import read_vtk_network
from collections import defaultdict
from itertools import repeat

def map_vasculature(fname):

    # Read network information from 'fname', and extract radii and
    # roots information as FEniCS cell and vertex functions,
    # respectively.
    mesh, radii, roots = read_vtk_network(fname)

    # Use networkx directed graph representation of mesh for convenience.
    G = nx.DiGraph()
    G.add_edges_from(mesh.cells())
    degrees = np.array([G.degree(i) for i in range(G.number_of_nodes())])
    print("Vertices where degree == 1: ", np.where(degrees == 1)[0])

    # Reduce the vascular network to a bifurcating tree and keep a map
    # from tree element index to vessel segment indices.

    # Start with the root node with lowest vertex index.
    first_roots = np.where(roots.array() == 2)[0]
    print("Vertices where root == 2: ", first_roots)
    more_roots = np.where(roots.array() == 1)[0]
    print("Vertices where root == 1: ", more_roots)
    #print(mesh.cells()[0])
        
    exit()

    # Reduce each subnetwork to bifurcating tree
    #num_edges = mesh.num_cells()
    #for c in dolfin.cells(mesh):
    #    print(c.index())
    #print("mesh.num_cells() = ", mesh.num_cells())
    
    print(mesh.cells())
    print(radii.array())
    print(roots.array())

    #print(roots, num_roots)

    
if __name__ == '__main__':

    map_vasculature("../mesh/networks/arteries_smooth.vtk")

    exit()
