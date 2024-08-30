from fenics import *
from xii import *
import numpy_indexed as npi
import numpy as np
from IPython import embed

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

def map_on_global(uh, parentmesh, eps_digits=12):
    # map the solution back on the whole domain
    el = uh.function_space().ufl_element()
    V = uh.function_space()
    if el.family()=="Lagrange":
        bdim = 1 if len(uh.ufl_shape)==0 else uh.ufl_shape[0]
    if el.family()=="Brezzi-Douglas-Marini":
        bdim = 1
    if el.family()=='Discontinuous Lagrange' and el.degree()==0:
        bdim = 1
    V_glob   = FunctionSpace(parentmesh, el)
    c_coords = np.round(V.tabulate_dof_coordinates()[::bdim,:], eps_digits)
    p_coords = np.round(V_glob.tabulate_dof_coordinates()[::bdim,:], eps_digits)
    idxmap = npi.indices(p_coords, c_coords, axis=0)
    uh_global = Function(V_glob)
    for i in range(bdim):
        uh_global.vector()[idxmap*bdim + i] = uh.vector()[i::bdim]
    assert np.isclose(uh.vector().sum(), uh_global.vector().sum())
    return uh_global

"""

mesh = Mesh()
with XDMFFile(f'mesh/T1/volmesh/mesh.xdmf') as f:
    f.read(mesh)
    gdim = mesh.geometric_dimension()
    label = MeshFunction('size_t', mesh, gdim, 0)
    f.read(label, 'label')

emb = EmbeddedMesh(label, [CSFID,LVID,V34ID]) 
"""
mesh = UnitSquareMesh(4,4)
left = CompiledSubDomain("x[0] <= 0.5")
mf = MeshFunction("size_t", mesh, 2, 0)
left.mark(mf, 1)
emb = EmbeddedMesh(mf, 1) 

gdim = mesh.geometric_dimension()
CG2emb = FunctionSpace(emb, "BDM", 1)
uh = interpolate(Expression(("sin(x[0]*pi)", "sin(x[1]*pi)"),degree=1), CG2emb)
bc = DirichletBC(CG2emb, Constant([0]*gdim), "on_boundary")
#bc.apply(uh.vector())

uh_global = map_on_global(uh, mesh)

print(assemble(div(uh)*div(uh)*dx))
print(assemble(div(uh_global)*div(uh_global)*dx))

uhdg = interpolate(uh, VectorFunctionSpace(emb, "DG", 1))
with XDMFFile("uh.xdmf") as file:
    file.write_checkpoint(uhdg, "uh")

uh.set_allow_extrapolation(True)
uh_global_dg = interpolate(uh, VectorFunctionSpace(mesh, "DG", 1))
with XDMFFile("uh_gobal.xdmf") as file:
    file.write_checkpoint(uh_global_dg, "uh_global")
