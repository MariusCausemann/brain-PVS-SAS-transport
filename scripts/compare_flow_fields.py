from fenics import *

def read_velocity(vel_file, vel_mesh=None):
    with HDF5File(MPI.comm_world, vel_file,'r') as f:
        if vel_mesh is None:
            vel_mesh = Mesh()
            f.read(vel_mesh, "mesh", False)
        v_elem = eval(f.attributes("/velocity").to_dict()["signature"])
        V = FunctionSpace(vel_mesh, v_elem)
        velocity_csf = Function(V)
        f.read(velocity_csf, "velocity")
    return velocity_csf

vel_file1 = "results/csf_flow/sas_flow/flow.hdf"
vel_file2 = "results/csf_flow/sas_flow_pointsource/flow.hdf"

u1 = read_velocity(vel_file1)
u2 = read_velocity(vel_file2, vel_mesh=u1.function_space().mesh())
diff = u1 - u2
L2_err = sqrt(assemble(inner(diff, diff)*dx))
L2_u1 = sqrt(assemble(inner(u1, u1)*dx))

L2_err_rel = L2_err / L2_u1

V = u1.function_space()
error_field = Function(V)

error_field.vector()[:] = u1.vector() - u2.vector()
File("velocity_error.pvd") << error_field