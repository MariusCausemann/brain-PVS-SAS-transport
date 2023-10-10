from petsc4py import PETSc
from dolfin import *
# import xii


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

    vein = Mesh('../mesh/venous_network.xml')
    vein_radii = MeshFunction('double', vein, '../mesh/venous_network_radii.xml')
    vein_radii = as_P0_function(vein_radii)

    file = File("../mesh/vein.pvd")

    file << vein_radii
    artery = Mesh('../mesh/arterial_network.xml')
    artery_radii = MeshFunction('double', artery, '../mesh/arterial_network_radii.xml')
    artery_radii = as_P0_function(artery_radii)

    file2 = File("../mesh/artery.pvd")

    file2 << artery_radii
#sas = Mesh()
#with XDMFFile('../mesh/volmesh/mesh.xdmf') as f:
#        f.read(sas)
