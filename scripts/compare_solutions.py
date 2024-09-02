from fenics import *

filefine = "results/csf_flow/sas_flow_HighRes/csf_v.xdmf"
filecoarse = "results/csf_flow/sas_flow/csf_v.xdmf"
space = "DG"

mesh1 = Mesh()
with XDMFFile(filefine) as f:
    f.read(mesh1)
    V1 = VectorFunctionSpace(mesh1, space, 1)
    sol_fine = Function(V1)
    f.read_checkpoint(sol_fine, "velocity")

mesh2 = Mesh()
with XDMFFile(filecoarse) as f:
    f.read(mesh2)
    V2 = VectorFunctionSpace(mesh2, space, 1)
    sol_coarse = Function(V2)
    f.read_checkpoint(sol_coarse, "velocity")

sol_coarse.set_allow_extrapolation(True)
coarse_on_fine = project(sol_coarse, V1, solver_type="mumps")

diff = Function(V1)
diff.vector()[:] = sol_fine.vector()[:] - coarse_on_fine.vector()[:]

L2fine = norm(sol_fine)
L2coarse =  norm(coarse_on_fine)
L2diff = norm(diff)

relerror = L2diff / L2fine