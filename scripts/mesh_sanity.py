import dolfin as df
import networkx as nx
import numpy as np


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
    opts.setValue('eps_monitor', None)
    opts.setValue('eps_view', None)
    opts.setValue('eps_view_pre', None)    
    opts.setValue('eps_tol', 1E-8)
    opts.setValue('eps_type', 'krylovschur')
    opts.setValue('st_ksp_rtol', 1E-12)
    opts.setValue('st_ksp_type', 'cg')
    opts.setValue('st_pc_type', 'hypre')
    opts.setValue('st_ksp_monitor_true_residual', None)
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
        pairs.append((val, foo_mode))
    pairs = sorted(pairs, key=lambda p: p[0])
