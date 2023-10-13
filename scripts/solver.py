from petsc4py import PETSc
from dolfin import *
import xii  
import pyvista as pv


def read_vtk_network(filename):

    netw = pv.read(filename)
    mesh = Mesh()
    ed = MeshEditor()
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
    radii = MeshFunction('double', mesh, 1, 0)
    roots = MeshFunction('size_t', mesh, 0, 0)

    roots.array()[:] = netw["root"]
    netw = netw.point_data_to_cell_data()
    radii.array()[:] = netw["radius"]

    return mesh, radii, roots


def pcws_constant(subdomains, values):
    mesh = subdomains.mesh()
    dx = Measure('dx', domain=mesh, subdomain_data=subdomains)

    V = FunctionSpace(mesh, 'DG', 0)
    v = TestFunction(V)
    hK = CellVolume(mesh)

    form = sum((1/hK)*inner(Constant(value), v)*dx(tag) for tag, value in values.items())

    foo = Function(V)
    assemble(form, foo.vector())

    return foo

def as_P0_function(mesh_f):
    '''Represent as DG0'''
    mesh = mesh_f.mesh()
    
    assert mesh_f.dim() == mesh.topology().dim()
    P0 = FunctionSpace(mesh, 'DG', 0)
    f = Function(P0)
    f.vector().set_local(mesh_f.array())

    return f


def ksp_mat(tensor):
    '''Underlying PETSc thing'''
    return as_backend_type(tensor).mat()


def ksp_vec(tensor):
    '''Underlying PETSc thing'''
    return as_backend_type(tensor).vec()

# --------------------------------------------------------------------


if __name__ == '__main__':

    vein, vein_radii, vein_roots = read_vtk_network("../mesh/networks/venes_smooth.vtk")
    vein_radii = as_P0_function(vein_radii)

    artery, artery_radii, artery_roots = read_vtk_network("../mesh/networks/arteries_smooth.vtk")
    artery_radii = as_P0_function(artery_radii)
    

    sas = Mesh()
    with XDMFFile('../mesh/volmesh/mesh.xdmf') as f:
        f.read(sas)

    V = FunctionSpace(sas, 'CG', 1)
    Qa = FunctionSpace(artery, 'CG', 1)
    Qv = FunctionSpace(vein, 'CG', 1)
    W = [V, Qa, Qv]

    u, pa, pv = map(TrialFunction, W)
    v, qa, qv = map(TestFunction, W)

    # Things for restriction
    dx_a = Measure('dx', domain=artery)
    artery_shape = xii.Circle(radius=artery_radii, degree=20)
    ua, va = (xii.Average(x, artery, artery_shape) for x in (u, v))

    dx_v = Measure('dx', domain=vein)
    vein_shape = xii.Circle(radius=vein_radii, degree=20)
    uv, vv = (xii.Average(x, vein, vein_shape) for x in (u, v))

    a = xii.block_form(W, 2)
    a[0][0] = inner(grad(u), grad(v))*dx + inner(u, v)*dx + inner(ua, va)*dx_a + inner(uv, vv)*dx_v
    a[0][1] = -inner(pa, va)*dx_a
    a[0][2] = -inner(pv, vv)*dx_v

    a[1][0] = -inner(qa, ua)*dx_a
    a[1][1] = inner(grad(pa), grad(qa))*dx + inner(pa, qa)*dx

    a[2][0] = -inner(qv, uv)*dx_v
    a[2][2] = inner(grad(pv), grad(qv))*dx + inner(pv, qv)*dx

    L = xii.block_form(W, 1)
    

    V_bcs  =  []
    Qa_bcs = [DirichletBC(Qa, Expression('x[0]+x[1]+x[2]', degree=1), 'on_boundary')]
    Qv_bcs = [DirichletBC(Qv, Expression('x[0]+x[1]+x[2]', degree=1), 'on_boundary')]
    W_bcs = [V_bcs, Qa_bcs, Qv_bcs]

    AA, bb = map(xii.ii_assemble, (a, L))
    A, b = xii.apply_bc(AA, bb, bcs=W_bcs)
    A_, b_ = (ksp_mat(xii.ii_convert(A)), ksp_vec(xii.ii_convert(b)))

    ksp = PETSc.KSP().create()

    opts = PETSc.Options()
    opts.setValue('ksp_type', 'cg')    
    opts.setValue('ksp_view', None)
    opts.setValue('ksp_view_eigenvalues', None)
    opts.setValue('ksp_converged_reason', None)
    opts.setValue('ksp_monitor_true_residual', None)
    opts.setValue('ksp_rtol', 1E-40)
    opts.setValue('ksp_atol', 1E-12)
    opts.setValue('pc_type', 'hypre')
    opts.setValue('ksp_initial_guess_nonzero', 1)

    ksp.setOperators(A_, A_)
    ksp.setFromOptions()
    print('Start solve')
    wh = xii.ii_Function(W)
    
    x_ = b_.copy()
    x_ *= 0 
    ksp.solve(b_, x_)
    # NOTE: solve(b_, ksp_vec(wh.vector())) segfault most likely because
    # of type incompatibility seq is expected and we have nest
    wh.vector()[:] = PETScVector(x_)

    File('../results/uh_sas.pvd') << wh[0]
    File('../results/uh_artery.pvd') << wh[1]
    File('../results/uh_vein.pvd') << wh[2]

    import matplotlib.pyplot as plt
    from scipy.sparse import csr_matrix

    that = xii.ii_convert(AA[1][0])
    that = csr_matrix(ksp_mat(that).getValuesCSR()[::-1],
                      shape=(that.size(0), that.size(1)))
    
    fig, ax = plt.subplots()
    ax.spy(that)
    plt.show()
