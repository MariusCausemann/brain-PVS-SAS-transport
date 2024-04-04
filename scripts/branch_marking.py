import os.path
from collections import defaultdict
from itertools import chain, combinations
import networkx as nx
import dolfin as df
import numpy as np


def imap(mapping):
    '''Invert dict mapping item to collection'''
    inverse = defaultdict(set)
    [inverse[item].add(key) for key in mapping for item in mapping[key]]

    return inverse


def is_loop(mesh):
    '''No bifurcations'''
    assert mesh.topology().dim() == 1
    
    tdim = mesh.topology().dim()
    _, f2c = mesh.init(tdim-1, tdim), mesh.topology()(tdim-1, tdim)
    
    return all(len(f2c(f)) == 2 for f in range(mesh.num_vertices()))


def color_branches(mesh):
    '''Start/end is a terminal'''
    assert mesh.topology().dim() == 1

    _, v2c = mesh.init(0, 1), mesh.topology()(0, 1)
    c2v = mesh.topology()(1, 0)

    terminals = {v: set(v2c(v)) for v in range(mesh.num_vertices()) if len(v2c(v)) != 2}
    
    cell_f = df.MeshFunction('size_t', mesh, 1, 0)
    if not terminals:
        cell_f.set_all(1)
        return cell_f, [], [1]

    def next_vertex(c, v, c2v=c2v):
        v0, v1 = c2v(c)
        return v1 if v == v0 else v0

    def next_cell(v, c, v2c=v2c):
        c0, c1 = v2c(v)
        return c1 if c == c0 else c0
    
    values = cell_f.array()
    branch_colors, loop_colors, color = [], [], 0

    exhausted = False
    while not exhausted:
        vertex = max(terminals, key=lambda v: terminals[v])
        vertex_cells = terminals[vertex]

        exhausted = len(vertex_cells) == 0
        # The idea is to walk from vertex following the cell
        while vertex_cells:
            link_cell = vertex_cells.pop()
            v0 = vertex

            branch = [link_cell]
            # v0 --
            while next_vertex(link_cell, v0) not in terminals:
                # -- v0 ==
                v0 = next_vertex(link_cell, v0)
                # Because we have not terminal, ==
                link_cell = next_cell(v0, link_cell)
                branch.append(link_cell)
            # Once we reached the terminal
            v0 = next_vertex(link_cell, v0)

            color += 1
            if v0 == vertex:
                loop_colors.append(color)
            else:
                branch_colors.append(color)
            values[branch] = color
            
            # Preclude leaving from vertex in a look
            link_cell in vertex_cells and vertex_cells.remove(link_cell)
            # If we arrived to some other terminal, we don't want to leave from it by the
            # same way we arrived
            v0 in terminals and link_cell in terminals[v0] and terminals[v0].remove(link_cell)

    return cell_f, branch_colors, loop_colors


def walk_vertices(arg, tag=None, is_loop=False):
    '''Walk vertices in a linked way'''
    assert isinstance(arg, df.Mesh) or isinstance(arg, df.cpp.mesh.MeshFunctionSizet)
    # Branch
    if isinstance(arg, df.cpp.mesh.MeshFunctionSizet):
        mesh = arg.mesh()
        assert arg.dim() == 1
    else:
        mesh = arg
        
    assert mesh.topology().dim() == 1 and mesh.geometry().dim() > 1
    # The boring cese
    assert mesh.num_cells() > 1
    
    c2v = mesh.topology()(1, 0)
    _, v2c = mesh.init(0, 1), mesh.topology()(0, 1)

    cells = walk_cells(arg.array(), tag=tag, c2v=c2v, v2c=v2c, is_loop=is_loop)
    cell, orient = next(cells)

    vertices = c2v(cell) if orient else reversed(c2v(cell))
    for v in vertices:
        yield v
    
    for cell, orient in cells:
        yield list(c2v(cell) if orient else reversed(c2v(cell)))[-1]
        

def walk_cells(cell_f, tag, c2v, v2c, is_loop):
    '''Walk cells where cell_f == tag in a linked way'''
    cell_indices, = np.where(cell_f == tag)
    # Localize to tags
    c2v = {c: c2v(c) for c in cell_indices}
    v2c = imap(c2v)

    # We return cell index together with orientation, i.e. True if link
    # is v0, v1 False if link is v1, v0
    def next_vertex(c, v, c2v=c2v):
        v0, v1 = c2v[c]
        return v1 if v == v0 else v0

    def next_cell(v, c, v2c=v2c):
        c0, c1 = v2c[v]
        return c1 if c == c0 else c0

    if is_loop:
        # Pick first marked cell
        link_cell = cell_indices[0]
        # For loop we pick where to start as either of the first cell
        start, v1 = c2v[link_cell]
        # ... and we terminate once we reach the start again
        end = start
    else:
        # If this is a branch we need two end cells/vertices
        # One is a start the other is end
        start, end = [v for v in v2c if len(v2c[v]) == 1]
        
        link_cell, = v2c[start]
        # The linking vertex is not the start
        v1,  = set(c2v[link_cell]) - set((start, ))
        
    yield link_cell, c2v[link_cell][-1] == v1

    v0 = start
    while next_vertex(link_cell, v0) != end:
        # -- v0 ==
        v0 = next_vertex(link_cell, v0)
        # Because we have not terminal, ==
        link_cell = next_cell(v0, link_cell)

        yield link_cell, c2v[link_cell][0] == v0


def minimal_color(cell_f):
    '''A node is the branch colored by some color'''
    mesh = cell_f.mesh()
    assert mesh.topology().dim() == 1
    cell_f_values = cell_f.array()
    
    _, v2c = mesh.init(0, 1), mesh.topology()(0, 1)
    c2v = mesh.topology()(1, 0)

    branching_cells = {v: set(v2c(v)) for v in range(mesh.num_vertices()) if len(v2c(v)) > 2}

    graph_edges = set(
        chain(*(combinations([cell_f_values[c] for c in cs], 2) for cs in branching_cells.values()))
    )
    
    g = nx.Graph()
    g.add_edges_from(graph_edges)

    new_coloring = nx.algorithms.greedy_color(g)

    new_cell_f = df.MeshFunction('size_t', mesh, 1, 0)
    new_cell_f_values = new_cell_f.array()

    for old, new in new_coloring.items():
        new_cell_f_values[cell_f_values == old] = new
    return new_cell_f


def first(iterable):
    '''The[0]'''
    return next(iter(iterable))

# --------------------------------------------------------------------

if __name__ == '__main__':
    from solver import read_vtk_network
    from collections import defaultdict
    from itertools import repeat

    # Place output in separate subdirectory:
    subdir = "branches"
    
    artery, artery_radii, artery_roots = read_vtk_network("../mesh/networks/arteries_smooth.vtk", rescale_mm2m=False)

    # Find branch and mark the cells that make it up by unique color
    marking_branch, branch_colors, loop_colors = color_branches(artery)
    # This may be too many colors. An alternative is to minimally color the graph
    df.File(os.path.join(subdir, 'unique_branch.pvd')) << marking_branch

    mesh = marking_branch.mesh()
    Q = df.FunctionSpace(mesh, 'DG', 0)
    q = df.TestFunction(Q)
    hK = df.CellVolume(mesh)

    foo = df.Function(Q)   # Holds in edges of branch the lenght of the branch
    foo_values = foo.vector().get_local()
    
    active_cell_f = df.MeshFunction('size_t', mesh, 1, 0)
    dx_active = df.Measure('dx', domain=mesh, subdomain_data=active_cell_f)
    dx_marked = df.Measure('dx', domain=mesh, subdomain_data=marking_branch)
    
    length_form = q*dx_active(1)

    sanity_check = False  # Compare manual on/off marking with the original coloring
    
    active_cell_f_array = active_cell_f.array()
    cell_colors = marking_branch.array()
    for color in branch_colors:
        active_cell_f_array[cell_colors == color] = 1

        branch_lengths = df.assemble(length_form)
        branch_length = branch_lengths.sum() 
        foo_values[branch_lengths.get_local() > 0] = branch_length
        # NOTE: this will trigget the compiler a lot, that's why we do
        # logic above. But I keep it here for reference anyways
        if sanity_check:
            target = df.assemble(df.Constant(1)*dx_marked(color))
            assert target > 0 and abs(branch_length - target) < 1E-13, ((color, target, branch_length))
        
        active_cell_f_array *= 0
    foo.vector().set_local(foo_values)
    with df.XDMFFile(os.path.join(subdir, 'branch_length_read.xdmf')) as xdmf:
        xdmf.write_checkpoint(foo, "branch_length")
    df.File(os.path.join(subdir, 'branch_length.pvd')) << foo
    
    minimal_marking_branch = minimal_color(marking_branch)
    df.File(os.path.join(subdir, 'minimal_branch.pvd')) << minimal_marking_branch
    
    # Getting the network out in different formats
    mesh = marking_branch.mesh()
    x = mesh.coordinates()

    _, c2v = mesh.init(1, 0), mesh.topology()(1, 0)
    _, v2c = mesh.init(0, 1), mesh.topology()(0, 1)
    
    g = nx.Graph()
    g.add_edges_from(mesh.cells())
    # Topology of branches; vertex -> colors
    color_connectivity = defaultdict(list)   # C
    branch_paths, branch_radii = {}, {}
        
    radii = artery_radii.array()
    edge_colors = marking_branch.array()
    # Since each branch is colored we can visit it as part from source to destination
    for color in branch_colors:
        marked_edges, = np.where(edge_colors == color)
        # Encode vertex to edge of the branch to find terminals
        branch_connectivity = defaultdict(list)
        # Start and stop
        for edge in marked_edges:
            v0, v1 = c2v(edge)
            branch_connectivity[v0].append(edge)
            branch_connectivity[v1].append(edge)            
        source, dest = [node for (node, cells) in branch_connectivity.items() if len(cells) == 1]

        color_connectivity[source].append(color)
        color_connectivity[dest].append(color)
        
        path = nx.algorithms.shortest_path(g, source, dest)
        # Here the path is in terms of nodes. We represent it as edges ...
        branch_paths[color] = tuple(zip(path[:-1], path[1:]))
        # and for each edge we grab the radius
        branch_radii[color] = tuple(radii[first(set(v2c(v0)) & set(v2c(v1)))] for (v0, v1) in branch_paths[color])

    # Reduced graph is one where we discregard the interior geometry of the branch
    reduced_edges = defaultdict(list)
    for vertex, colors in color_connectivity.items():
        [reduced_edges[color].append(vertex) for color in colors]

    reduced_mesh = df.Mesh()
    editor = df.MeshEditor()
    editor.open(reduced_mesh, 'interval', 1, 3)
    editor.init_vertices(len(color_connectivity))
    editor.init_cells(len(reduced_edges))

    reorder = {}
    for new, old in enumerate(color_connectivity):
        editor.add_vertex(new, x[old])
        reorder[old] = new

    for idx, color in enumerate(reduced_edges):
        cell = tuple(reorder[v] for v in reduced_edges[color])
        editor.add_cell(idx, cell)
    editor.close()

    reduced_color_f = df.MeshFunction('size_t', reduced_mesh, 1, 0)
    reduced_color_f.array()[:] = np.fromiter(reduced_edges, dtype='uintp')

    df.File(os.path.join(subdir, 'reduced_branch.pvd')) << reduced_color_f
