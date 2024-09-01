from fenics import *

file1 = "results/csf_flow/sas_flow/csf_v.xdmf"
space = "DG"

mesh1 = Mesh()
with XDMFFile(file1) as f:
    f.read(mesh1)
    V = VectorFunctionSpace(mesh1, space, 1)
    velocity_sas = Function(V)
    f.read_checkpoint(velocity_sas, "velocity")