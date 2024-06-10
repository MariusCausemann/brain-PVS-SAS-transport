import dolfin as df
import networkx as nx
import numpy as np

from xii import *
import pyvista as pv 

def color_connected_components(mesh):
    '''
    Cell function with colors corresponding to connected components of mesh graph. 
    Largest components first.
    '''
    print(f'Coloring components in mesh with {mesh.num_cells()} cells')
    tdim = mesh.topology().dim()
    mesh.init(tdim-1, tdim)  #  Facet to cell

    g = nx.Graph()
    g.add_edges_from(tuple(facet.entities(tdim))
                     for facet in df.facets(mesh) if len(facet.entities(tdim)) > 1)

    ccs = sorted(list(nx.algorithms.connected_components(g)))

    
    components = df.MeshFunction('size_t', mesh, tdim, 0)
    values = components.array()
    for (color, component) in enumerate(ccs, 1):
        print(f'Component {color} with {len(component)} cells')
        values[list(component)] = color
    return components

# ----------------------------------------------------------------------

if __name__ == '__main__':
    from slepc4py import SLEPc
    from petsc4py import PETSc
    
<<<<<<< HEAD
    path = '/home/mirok/Downloads/mesh.xdmf'
    mesh = df.Mesh()
    with df.XDMFFile(mesh.mpi_comm(), path) as f:
        f.read(mesh)
    gdim = mesh.geometry().dim()

    # Extract mesh as largest connected component
    if mesh.mpi_comm().size == 1:
        cell_f = color_connected_components(mesh)        
        ncomps = len(np.unique(cell_f.array()))
        if ncomps > 1:
            mesh = df.SubMesh(mesh, cell_f, 1)
=======
    path = 'mesh/mid_mesh/volmesh/mesh.xdmf' 
 
    # get mesh 
    sas = df.Mesh()
    with df.XDMFFile('mesh/mid_mesh/volmesh/mesh.xdmf') as f:
        f.read(sas)
    gdim = sas.geometric_dimension()
    sas_components = df.MeshFunction('size_t', sas, gdim, 0)
    f.read(sas_components, 'sas_components')
    sas.scale(1e-3)  # scale from mm to m

    # obtain the sas mesh 
    sas_outer = EmbeddedMesh(sas_components, 1) 
    # create boundary markers 
    boundary_markers = df.MeshFunction("size_t", sas, sas.topology().dim() - 1)
    boundary_markers.set_all(0)
    # Sub domain for efflux route (mark whole boundary of the full domain) 
    class Efflux(df.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary # and (pow(x[0] - 0.16122659231266656,2) +pow(x[1] - 0.18953817553340035,2) + pow(x[2] - 0.1823329158961796 ,2 ) > pow(.09,2))
    efflux = Efflux()
    efflux.mark(boundary_markers, 1)

    # translate markers to the sas outer mesh 
    boundary_markers_outer = df.MeshFunction("size_t", sas_outer, sas_outer.topology().dim() - 1, 0)
    df.DomainBoundary().mark(boundary_markers_outer, 2)
    sas_outer.translate_markers(boundary_markers, (1, ), marker_f=boundary_markers_outer)

    #File("fpp.pvd") << boundary_markers_outer


    # Extract mesh as largest connected component
    if sas_outer.mpi_comm().size == 1:
        cell_f = color_connected_components(sas_outer)        
        ncomps = len(np.unique(cell_f.array()))
        if ncomps > 1:
            mesh = df.SubMesh(sas_outer, cell_f, 1)

    df.XDMFFile("mesh/mid_mesh/volmesh/sas_outer.xdmf").write(cell_f) 

    sas_outer = EmbeddedMesh(cell_f, 1) 
>>>>>>> f4126d1a02f5c4a0a7888520ba351c54b0b10885

    df.set_log_level(10)
    # Sanity check in terms of vector Poisson. Don't want any zero
    # eigenvalues
    V = df.VectorFunctionSpace(mesh, 'CG', 1)
    print(f'dim(V) = {V.dim()}')
    bcs = df.DirichletBC(V, df.Constant((0, )*gdim), 'on_boundary')
    
    u, v = df.TrialFunction(V), df.TestFunction(V)

    a = df.inner(df.grad(u), df.grad(v))*df.dx
    m = df.inner(u, v)*df.dx
    L = df.inner(df.Constant((0, )*gdim), v)*df.dx

    A, _ = df.assemble_system(a, L, bcs)
    B, _ = df.assemble_system(m, L, bcs)

    A, B = (df.as_backend_type(mat).mat() for mat in (A, B))

    opts = PETSc.Options()
    dopts = opts.getAll()
    
    opts.setValue('eps_max_it', 50_000)
    opts.setValue('eps_nev', 10)
<<<<<<< HEAD
    opts.setValue('eps_monitor', None)
=======
    #opts.setValue('eps_monitor', None)
>>>>>>> f4126d1a02f5c4a0a7888520ba351c54b0b10885
    opts.setValue('eps_view', None)
    opts.setValue('eps_view_pre', None)    
    opts.setValue('eps_tol', 1E-8)
    opts.setValue('eps_type', 'krylovschur')
    opts.setValue('st_ksp_rtol', 1E-12)
    opts.setValue('st_ksp_type', 'cg')
    opts.setValue('st_pc_type', 'hypre')
<<<<<<< HEAD
    opts.setValue('st_ksp_monitor_true_residual', None)

    #opts.setValue('st_ksp_monitor_true_residual', None)
    # Setup the eigensolver
    E = SLEPc.EPS().create()
    E.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    E.setOperators(A, B)
    E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)

    E.setFromOptions()
    E.solve()

    its = E.getIterationNumber()
    nconv = E.getConverged()

    pairs = []
    for i in range(nconv):
        mode = A.createVecLeft()
        val = E.getEigenpair(i, mode)
        print(f'Eigenvalue {i} is {val}')
        pairs.append((val, mode))
    pairs = sorted(pairs, key=lambda p: p[0])
