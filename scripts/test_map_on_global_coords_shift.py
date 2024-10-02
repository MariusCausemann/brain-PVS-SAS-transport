from fenics import *
from xii import *
import numpy as np
import numpy_indexed as npi
from pykdtree.kdtree import KDTree

def map_kdtree(data_points, query_points):
    tree = KDTree(data_points)
    dist, idx = tree.query(query_points, k=1)
    return idx

def get_global_space(uh, mesh):
    el = uh.function_space().ufl_element()
    return FunctionSpace(mesh, el)

def map_dg_on_global(uh, uh_global=None, parentmesh=None):
    # map the solution back on the whole domain
    childmesh = uh.function_space().mesh()
    if uh_global is not None:
        parentmesh = uh_global.function_space().mesh()
    elif parentmesh is not None:
        uh_global = Function(get_global_space(uh, parentmesh))
    eldim = 1 if len(uh.ufl_shape)==0 else uh.ufl_shape[0]
    C = uh.function_space().tabulate_dof_coordinates()[::eldim,:]
    P = uh_global.function_space().tabulate_dof_coordinates()[::eldim,:]
    eldim = 1 if len(uh.ufl_shape)==0 else uh.ufl_shape[0]

    dofs_per_cell = len(uh.function_space().dofmap().cell_dofs(0))
    Cc = FunctionSpace(childmesh, "DG", 0).tabulate_dof_coordinates()
    Pc = FunctionSpace(parentmesh, "DG", 0).tabulate_dof_coordinates()
    # shift coords towards the cell center
    C_shift = 0.5*(C + np.repeat(Cc, int(dofs_per_cell/eldim), axis=0))
    P_shift = 0.5*(P + np.repeat(Pc, int(dofs_per_cell/eldim), axis=0))

    #idxmap = npi.indices(np.round(P_shift, ndigits), np.round(C_shift, ndigits), axis=0)
    idxmap = map_kdtree(P_shift, C_shift)
    for i in range(eldim):
        uh_global.vector()[idxmap*eldim + i] = uh.vector()[i::eldim]
    return uh_global
CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

if __name__=="__main__":
    mesh = UnitCubeMesh(20,20,20)
    gdim = mesh.geometric_dimension()
    mf = MeshFunction("size_t", mesh, gdim, 0)
    left = CompiledSubDomain("x[0] <= 0.5")
    left.mark(mf, 1)
    emb = EmbeddedMesh(mf, [1])

    f = Expression(["x[0]"]*gdim, degree=2)

    V = VectorFunctionSpace(mesh, "DG", 1)
    Vemb  = VectorFunctionSpace(emb, "DG", 1)

    uh = interpolate(f, Vemb)
    uh_global = Function(V)

    uh_global = map_dg_on_global(uh, uh_global)
    uh_global = map_dg_on_global(uh, parentmesh=mesh)

    with XDMFFile("uh.xdmf") as file:
        file.write_checkpoint(uh, "uh")

    with XDMFFile("uh_global.xdmf") as file:
        file.write_checkpoint(uh_global, "uh")

    print(assemble(inner(f - uh, f-uh)*dx))

    dX = Measure('dx', domain=mesh, subdomain_data=mf)
    print(assemble(inner(f - uh_global, f-uh_global)*dX(1)))


