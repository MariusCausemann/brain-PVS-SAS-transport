from fenics import *
from xii import *
import numpy as np
import numpy_indexed as npi

def map_on_global(uh, uh_global):
    # map the solution back on the whole domain
    childmesh = uh.function_space().mesh()
    parentmesh = uh_global.function_space().mesh()
    cell_map = childmesh.parent_entity_map[parentmesh.id()][childmesh.geometric_dimension()]
    C = uh.function_space().tabulate_dof_coordinates()
    P = uh_global.function_space().tabulate_dof_coordinates()
    for child, parent in cell_map.items():
        child_dofs = uh.function_space().dofmap().cell_dofs(child)
        parent_dofs = uh_global.function_space().dofmap().cell_dofs(parent)
        perm_map = npi.indices(C[child_dofs],P[parent_dofs], axis=0)
        uh_global.vector().vec().array_w[parent_dofs] = uh.vector().vec().array_r[child_dofs][perm_map]
        assert np.isclose(np.linalg.norm(C[child_dofs[perm_map]] - P[parent_dofs]) , 0)


def map_on_global_vec(uh, parentmesh):
    V = uh.function_space()
    Vsubparent = FunctionSpace(mesh, V.sub(0).ufl_element())
    uiglobs = [Function(Vsubparent) for v in V.split()]

    for ui, uigl in zip(uh.split(deepcopy=True), uiglobs):
        map_on_global(ui, uigl)
    return as_vector(uiglobs)


mesh = UnitCubeMesh(20,20,20)
gdim = mesh.geometric_dimension()
mf = MeshFunction("size_t", mesh, gdim, 0)
left = CompiledSubDomain("x[0] <= 0.5")
left.mark(mf, 1)
emb = EmbeddedMesh(mf, [1])

f = Expression(["x[0]"]*gdim, degree=2)

uh = project(f, VectorFunctionSpace(emb, "DG", 1))

uh_global = map_on_global_vec(uh, mesh)

with XDMFFile("uh.xdmf") as file:
    file.write_checkpoint(uh, "uh")

V = VectorFunctionSpace(mesh, "DG", 1)
Vemb  = VectorFunctionSpace(emb, "DG", 1)
with XDMFFile("uh_global.xdmf") as file:
    file.write_checkpoint(project(uh_global, V), "uh")

print(assemble(inner(f - uh, f-uh)*dx))

dX = Measure('dx', domain=mesh, subdomain_data=mf)
print(assemble(inner(f - uh_global, f-uh_global)*dX(1)))


