from fenics import *
from xii import *
import numpy_indexed as npi
import numpy as np

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

def map_on_global(uh, parentmesh, eps_digits=10):
    # map the solution back on the whole domain
    el = uh.function_space().ufl_element()
    V = uh.function_space()
    eldim = 1 if len(uh.ufl_shape)==0 else uh.ufl_shape[0]
    V_glob = FunctionSpace(parentmesh, el)
    c_coords = np.round(V.tabulate_dof_coordinates()[::eldim,:], eps_digits)
    p_coords = np.round(V_glob.tabulate_dof_coordinates()[::eldim,:], eps_digits)
    idxmap = npi.indices(p_coords, c_coords, axis=0)
    uh_global = Function(V_glob)
    for i in range(eldim):
        uh_global.vector()[idxmap*eldim + i] = uh.vector()[i::eldim]
    np.isclose(uh.vector().sum(), uh_global.vector().sum())
    return uh_global

mesh = Mesh()
with XDMFFile(f'mesh/T1/volmesh/mesh.xdmf') as f:
    f.read(mesh)
    gdim = mesh.geometric_dimension()
    label = MeshFunction('size_t', mesh, gdim, 0)
    f.read(label, 'label')

emb = EmbeddedMesh(label, [CSFID,LVID,V34ID]) 

CG2emb = VectorFunctionSpace(emb, "CG", 2)
uh = interpolate(Expression(("1", "0", "0"), degree=2), CG2emb)
bc = DirichletBC(CG2emb, Constant([0,0,0]), "on_boundary")
bc.apply(uh.vector())

uh_global = map_on_global(uh, mesh)

print(assemble(div(uh)*div(uh)*dx))
print(assemble(div(uh_global)*div(uh_global)*dx))

with XDMFFile("uh.xdmf") as file:
    file.write_checkpoint(uh, "uh")

with XDMFFile("uh_gobal.xdmf") as file:
    file.write_checkpoint(uh_global, "uh_global")
