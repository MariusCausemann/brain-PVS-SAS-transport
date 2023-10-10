
### rami: works with book-env: need to fix the environments mismatch !! 


import pyvista as pv
import numpy as np
from pyvtk import UnstructuredGrid, PointData, Scalars, VtkData
from collections import namedtuple
import networkx as nx
import numpy as np
import itertools, os
from vtk import *
import dolfin as df 
import vtk 
from vtk.util.numpy_support import vtk_to_numpy

SWCNeuron = namedtuple('swc', ('file', 'points', 'radii', 'graph'))


def read_swc(swc_file):
    '''SWC file to SWCNeuron'''
    points, radii, indices = [], [], [1]
    # 
    parse = lambda l: [fi(li) for fi, li in zip((int, int, float, float, float, float, int), l)]
    # 1 1 0.0 0.0 0.0 7.3875 -1
    with open(swc_file, 'r') as swc:
        # Get comments out of the way
        iter_points = map(lambda l: parse(l.strip().split()),
                          itertools.dropwhile(lambda l: l.startswith('#') or len(l) < 2, swc))
        # This is the first point - can't make pairs yet
        index, ntype, x, y, z, radius, parent = next(iter_points)
        # SWC starts numbering at 1
        assert index == 1, index

        points.append((x, y, z))
        radii.append(radius)

        G = nx.Graph()
        for line in iter_points:
            if not line:
                continue
            
            (index, ntype, x, y, z, radius, parent) = line

            points.append((x, y, z))
            radii.append(radius)
            indices.append(index)
            
            G.add_edge(index, parent)
    # A sanity check of the input is that the neuron is one graph
    assert len(list(nx.algorithms.connected_components(G))) == 1
    print('Neuron defined by %d points and %d edges' % (G.number_of_nodes(),
                                                        G.number_of_edges()))

    indices = np.array(indices) - 1

    points_ = np.array(points, dtype=float)
    points = np.zeros_like(points_)
    points[indices] = points_
    
    radii_ = np.array(radii, dtype=float)
    radii = np.zeros_like(radii_)
    radii[indices] = radii_
        
    return SWCNeuron(swc_file, points, radii, G)


def read_vtp(path):
    '''Dictionary with data in vtk file'''
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()

    data = reader.GetOutput().GetPointData()
    field_count = data.GetNumberOfArrays()
    data = {data.GetArrayName(i): vtk_to_numpy(data.GetArray(i)) for i in range(field_count)}
    
    if not data:
        data = reader.GetOutput().GetCellData()
        field_count = data.GetNumberOfArrays()
        data = {data.GetArrayName(i): vtk_to_numpy(data.GetArray(i)) for i in range(field_count)}
        
    # Add coordinates
    xyz = reader.GetOutput().GetPoints()
    points = vtk_to_numpy(xyz.GetData())

    # Cells
    output = reader.GetOutput()
    ncells = output.GetNumberOfCells()
    assert all(output.GetCellType(cell) == 3 for cell in range(ncells))
    
    cells = []
    for lid in range(ncells):
        cell = output.GetCell(lid)
        cells.append((cell.GetPointId(0), cell.GetPointId(1)))

    data['coordinates'] = points
    data['cells'] = np.array(cells)

    return data

def store_vtp(swc_neuron, vtp_file=''):
    # Store points in vtkPoints
    if not vtp_file:
        root, _ = os.path.splitext(swc_neuron.file)
        vtp_file = '.'.join([root, 'vtp'])
        
    coords = swc_neuron.points
    points = vtkPoints()
    for c in coords:
        points.InsertNextPoint(list(c))

    # Store edges in cell array
    lines = vtkCellArray()

    cell_data = []  # Midpoint radius
    for (v0, v1) in swc_neuron.graph.edges:
        line = vtkLine()
        line.GetPointIds().SetId(0, v0-1)
        line.GetPointIds().SetId(1, v1-1)
        lines.InsertNextCell(line)

        cell_data.append(0.5*(swc_neuron.radii[v0-1] + swc_neuron.radii[v1-1]))

    # Create a polydata to store 1d mesh in
    linesPolyData = vtkPolyData()
    linesPolyData.SetPoints(points)
    linesPolyData.SetLines(lines)

    # Write data from associated function
    data = vtkDoubleArray()
    data.SetName('radius')    
    data.SetNumberOfComponents(1)
        
        # store value of function at each coordinates
    for v in cell_data:
        data.InsertNextTuple([v])
    linesPolyData.GetCellData().AddArray(data)

    # Write to file
    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(vtp_file)
    writer.SetInputData(linesPolyData)
    writer.Update()
    writer.Write()

    return vtp_file




def dolfin_convert(data, mesh_f_key):
    '''Mesh and MesFunction<double> representation of data'''
    # In vtp as given in the dataset each segment has it own points. We
    # want them to be shared
    old_coordinates = list(map(tuple, data['coordinates']))
    old_vertex_data = data[mesh_f_key]
    
    vertex_to_idx = {vtx: idx for idx, vtx in enumerate(set(old_coordinates))}

    radius_vertex_data = len(old_vertex_data) != len(data['cells'])
    
    mesh_cells, cell_data = [], []
    # Cells in the mesh will be defined with respect to new ordering
    for k, (v0, v1) in enumerate(data['cells']):
        new0 = vertex_to_idx[old_coordinates[v0]]
        new1 = vertex_to_idx[old_coordinates[v1]]        
        mesh_cells.append((new0, new1))
        # We have vertex data for fenics in downsample to cell by mean
        if radius_vertex_data:
            cell_data.append(0.5*(old_vertex_data[v0] + old_vertex_data[v1]))
        else:
            cell_data.append(old_vertex_data[k])

    gdim, = set(map(len, old_coordinates))
    mesh_coordinates = np.zeros((len(vertex_to_idx), gdim))
    for vertex in vertex_to_idx:
        mesh_coordinates[vertex_to_idx[vertex]] = vertex

    # Make mesh
    mesh = df.Mesh()
    ed = df.MeshEditor()
    ed.open(mesh, 'interval', 1, gdim)
    ed.init_vertices(len(mesh_coordinates))
    ed.init_cells(len(mesh_cells))

    for vid, v in enumerate(mesh_coordinates):
        ed.add_vertex(vid, v)
    for cid, c in enumerate(mesh_cells):
        ed.add_cell(cid, c)
    ed.close()

    
    
    # Cell Function
    mesh_f = df.MeshFunction('double', mesh, 1, 0)
    mesh_f.array()[:] = np.array(cell_data)

    return mesh, mesh_f




neuron = read_swc('../mesh/arteries.swc')

print(store_vtp(neuron, '../mesh/arteries_from_swc.vtp'))
which = 'arteries_from_swc'
data = read_vtp('../mesh/arteries_from_swc.vtp')

#sas = Mesh()
#with XDMFFile('../mesh/volmesh/mesh.xdmf') as f:
#        f.read(sas)

 
mesh, mesh_f = dolfin_convert(data, mesh_f_key='radius')
assert mesh.hmin() > 0

df.File(f'../mesh/{which}.pvd') << mesh_f
# NOTE: save separately the mesh and mesh_f as XML for loading

df.File("../mesh/arterial_network.xml") << mesh 
df.File("../mesh/arterial_network_radii.xml") << mesh_f 




neuron = read_swc('../mesh/veins.swc')

print(store_vtp(neuron, '../mesh/veins_from_swc.vtp'))
data = read_vtp('../mesh/veins_from_swc.vtp')

#sas = Mesh()
#with XDMFFile('../mesh/volmesh/mesh.xdmf') as f:
#        f.read(sas)

 
mesh, mesh_f = dolfin_convert(data, mesh_f_key='radius')
assert mesh.hmin() > 0
which = 'veins_from_swc'
df.File(f'../mesh/{which}.pvd') << mesh_f
# NOTE: save separately the mesh and mesh_f as XML for loading

df.File("../mesh/venous_network.xml") << mesh 
df.File("../mesh/venous_network_radii.xml") << mesh_f 







