from itertools import chain, combinations
import networkx as nx
import dolfin as df
import numpy as np
import pymetis


def mesh_2_nxgraph(mesh, edge_weight_f=None):
    '''
    Cells[nodes in the graph] are connected by facets. This connectivity 
    defines the edges. Edge weights can be assigned to an edge by function 
    which takes in the two cell indices.
    '''
    tdim = mesh.topology().dim()
    fdim = tdim - 1

    _, c2f = mesh.init(tdim, fdim), mesh.topology()(tdim, fdim)    
    _, f2c = mesh.init(fdim, tdim), mesh.topology()(fdim, tdim)

    edges = chain(*(combinations(f.entities(tdim), 2) for f in df.facets(mesh)
                    if len(f.entities(tdim)) > 1))

    g = nx.Graph()    
    if edge_weight_f is None:
        g.add_edges_from(edges)
    else:
        make_edge = lambda edge: edge + (edge_weight_f(*edge), )
        g.add_weighted_edges_from(map(make_edge, edges))
    return g
                     

def nxgraph_2_metis(g, weight=''):
    '''Convert to structures that pymetis uses'''
    xadj, adjncy, eweights = [0], [], []
    for (this_node, connected_nodes) in sorted(g.adjacency()):
        xadj.append(xadj[-1] + len(connected_nodes))
        for node in connected_nodes:
            adjncy.append(node)
            if weight:
                edge_properties = connected_nodes[node]
                eweights.append(edge_properties[weight])
                
    metis_graph = {'xadj': xadj, 'adjncy': adjncy}
    # Edge weigths are obtained using the default naming the netwokx uses
    if eweights:
        eweights.extend(eweights)
        metis_graph['eweights'] = eweights

    return metis_graph


def partition(mesh, nsubs, weighted=False):
    '''
    Color mesh by mesh partitioner in `nsubs` subdomains. Return cell function
    colored by subdomains. With `weighted=True` we are trying to minimize the 
    sum of edge weigths involved in the cut - so large edge weight means less 
    likely split of the pair of cells.
    '''
    nxgraph = mesh_2_nxgraph(mesh, edge_weight_f=lambda u, v: u+v)
    
    metis_graph = nxgraph_2_metis(nxgraph, weight='weight' if weighted else '')
    ncuts, coloring = pymetis.part_graph(nsubs, **metis_graph)

    coloring = np.array(coloring)
    ncolors = len(np.unique(coloring))

    tdim = mesh.topology().dim()
    cell_f = df.MeshFunction('size_t', mesh, tdim, 1)
    values = cell_f.array()
    for color in range(ncolors):
        values[coloring == color] = color + 1

    return cell_f

# --------------------------------------------------------------------

if __name__ == '__main__':
    
    mesh = df.UnitSquareMesh(32, 32)

    cell_f = partition(mesh, 5, weighted=False)
    df.File('results/unweighted.pvd') << cell_f

