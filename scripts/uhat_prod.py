from dolfin import *
from solver import * 
import os
import numpy as np
from xii import *
from petsc4py import PETSc
import sympy as sp
import typer

def compute_avg_flow(csf_flow_model:  str):

    sas = Mesh()
    with XDMFFile(f'results/csf_flow/{csf_flow_model}/csf_v.xdmf') as f:
        f.read(sas)
        V = VectorFunctionSpace(sas, "CG", 3)
        velocity_sas = Function(V)
        f.read_checkpoint(velocity_sas, "velocity")
    
    #cell = sas.ufl_cell()
    #VP1_elem = VectorElement('Lagrange', cell, 1)
    #VP1 = FunctionSpace(sas, VP1_elem)
    #velocity_sas = Function(VP1) 
    #velocity_sas = interpolate(velocity_sas_ , VP1)
    # read the arterial network and mesh 
    pvs_ratio_artery  = 2.0 
    artery, artery_radii, artery_roots = read_vtk_network("mesh/networks/arteries_smooth.vtk", rescale_mm2m=False)
    artery_radii = as_Pk_function(artery_radii, k=3)

    # artery_radii.set_allow_extrapolation(True)
    pvs_radii = Function(artery_radii.function_space())
    pvs_radii.vector().set_local(pvs_ratio_artery*artery_radii.vector().get_local())
    pvs_radii.set_allow_extrapolation(True)
    artery_radii.set_allow_extrapolation(True)

    # define velocity average using Fenics_ii average functionality 
    disk_pvs  = Disk(radius    = pvs_radii, degree =11)
    disk_artery  = Disk(radius = artery_radii, degree =11)
    velocity_sas_averaged = Average(velocity_sas, artery, disk_pvs) -  Average(velocity_sas, artery, disk_artery)

    Qa = FunctionSpace(artery, 'CG', 1)
    tau = TangentCurve(artery)

    #  project velocity_sas_averaged \cdot tangent vector onto the function space over the network 

    u = TrialFunction(Qa)
    v = TestFunction(Qa)  
    dx_a = Measure('dx', domain=artery)
    a  = inner(u,v)*dx_a  
    L  = inner(dot(velocity_sas_averaged,tau), v)*dx_a
    A,b = map(ii_assemble, (a,L)) 
    A,b = map(ii_convert, (A,b)) 
    pvs_flow = Function(Qa) 

    solve(A, pvs_flow.vector(), b) 

    pvs_flow_vec = project(pvs_flow*tau, VectorFunctionSpace(artery, "DG", 0))

    with XDMFFile(f'results/csf_flow/{csf_flow_model}/pvs_flow.xdmf') as xdmf:
        xdmf.write_checkpoint(pvs_flow_vec, "velocity")

    File(f'results/csf_flow/{csf_flow_model}/pvs_flow.pvd') << pvs_flow_vec

if __name__ == "__main__":
    typer.run(compute_avg_flow)
