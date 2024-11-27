from dolfin import *
import xii
from solver import read_vtk_network, as_P0_function, get_mesh
import numpy as np
import os
from pvs_flow import pvs_flow_system
from test_map_on_global_coords_shift import map_kdtree, map_dg_on_global
from generate_synthseg_mesh import CSFID, V34ID, PARID, LVID, CSFNOFLOWID
import typer
from plotting_utils import read_config

m2mm = 1e3
m2mum = 1e6
Pa2mPa = 1e3

def map_pressure_on_network(p, pnet):
    idx = map_kdtree(p.function_space().tabulate_dof_coordinates(), 
                     pnet.function_space().tabulate_dof_coordinates())
    pnet.vector()[:] = p.vector()[idx]


def get_blood_flow_orientation(tau, radius_f, root_marker):
    mesh = tau.function_space().mesh()
    a, L, W = pvs_flow_system(radius_f, tau, 2, f=Constant(0))

    bcout = DirichletBC(W.sub(1), Constant(0), "on_boundary")
    bcin = DirichletBC(W.sub(1), Constant(1), root_marker, 2)

    A, b = assemble_system(a, L, [bcout, bcin])
    wh = Function(W)
    solve(A, wh.vector(), b)
    uh_mag, ph = wh.split(deepcopy=True)
    uh_sign = Function(uh_mag.function_space())
    uh_sign.vector()[:] = np.where(uh_mag.vector()[:] >=0, 1, -1)
    return uh_sign
# --------------------------------------------------------------------

def compute_pvs_flow(csf_flow_model:  str, network: str):

    radius_ratio = 2
    networkfile = f"mesh/networks/{network}_smooth.vtk"
    mesh, artery_radii, artery_roots = read_vtk_network(networkfile, rescale_mm2m=False)
    pressure_file = f"results/csf_flow/{csf_flow_model}/flow.hdf"
    results_dir = f"results/pvs_flow_prod/{csf_flow_model}-{network}/"
    csf_flow_config = read_config(f"configfiles/{csf_flow_model}.yml")
    radius_f = as_P0_function(artery_radii)
    
    pmesh = Mesh()
    with HDF5File(MPI.comm_world, pressure_file,'r') as f:
        f.read(pmesh, "mesh", False)
        print(f.attributes("/pressure").to_dict()["signature"])
        p_elem = eval(f.attributes("/pressure").to_dict()["signature"])
        DG = FunctionSpace(pmesh, p_elem)
        p = Function(DG)
        f.read(p, "pressure")

    _, sm = get_mesh(csf_flow_config["mesh"].split("_outer")[0] + ".xdmf")
    sas = xii.EmbeddedMesh(sm, [CSFID]) 
    par = xii.EmbeddedMesh(sm, [PARID, LVID, V34ID, CSFNOFLOWID]) 
    p = interpolate(p, FunctionSpace(sas, "DG", 1))

    p.set_allow_extrapolation(True)
    V = FunctionSpace(par, "CG", 1)
    u, v = TrialFunction(V), TestFunction(V)
    sol = Function(V)
    p_dirichletbc = DirichletBC(V, p, "on_boundary")
    solve(inner(grad(u), grad(v))*dx == Constant(0)*v*dx, sol, bcs=p_dirichletbc,
        solver_parameters={"linear_solver":"cg", "preconditioner":"amg"})

    p_combined = map_dg_on_global(p, parentmesh=sm.mesh())
    p_par = map_dg_on_global(interpolate(sol, FunctionSpace(par, "DG", 1)),
                             parentmesh=sm.mesh())

    p_combined.vector()[:] += p_par.vector()[:]

    p_combined.rename("p","p")


    # Grab the tangent of xii; bottom line is that this is vector valued
    # function on the network describing tangent to the edge. Orientation
    # is arbitrary as long same tau is used throught the code
    tau = xii.TangentCurve(mesh)
    downstream = get_blood_flow_orientation(tau, radius_f, artery_roots)

    a, L, W = pvs_flow_system(radius_f, downstream*tau, radius_ratio, f=Constant(0))

    pbc = Function(W.sub(1).collapse())
    map_pressure_on_network(p_combined, pbc)
    pbc.rename("p","p")

    bc = DirichletBC(W.sub(1), pbc, "on_boundary")

    A, b = assemble_system(a, L, [bc])

    wh = Function(W)
    solve(A, wh.vector(), b)

    uh_mag, ph = wh.split(deepcopy=True)

    pvs_flow_vec = project(uh_mag*tau*downstream, 
        VectorFunctionSpace(mesh, "DG", 1))
    
    ph.rename("p","p")
    pvs_flow_vec.rename("uh", "uh")
    os.makedirs(results_dir, exist_ok=True)

    with HDF5File(MPI.comm_world, f'{results_dir}/pvs_flow.hdf','w') as f:
        f.write(mesh, "mesh")
        f.write(pvs_flow_vec, "u")
        f.write(uh_mag, "umag")
        f.write(ph, "p")
        f.write(radius_f, "radii")

    with XDMFFile(f'{results_dir}/pvs_flow_vis.xdmf') as xdmf:
        xdmf.write(pvs_flow_vec, t=0)
        xdmf.write(ph, t=0)
    
if __name__ == "__main__":
    typer.run(compute_pvs_flow)
