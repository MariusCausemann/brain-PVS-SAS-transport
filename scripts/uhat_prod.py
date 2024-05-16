from dolfin import *
from solver import * 
import os
import numpy as np
from xii import *

# get mesh 
sas = Mesh()
with XDMFFile('mesh/T1/volmesh/mesh.xdmf') as f:
    f.read(sas)
    gdim = sas.geometric_dimension()
    vol_subdomains = MeshFunction('size_t', sas, gdim, 0)
    f.read(vol_subdomains, 'label')
    sas.scale(1e-3)  # scale from mm to m

# pick the sas mesh 
#sas_outer = EmbeddedMesh(vol_subdomains, 1) 

# Define function spaces for velocity and pressure
cell = sas.ufl_cell()  
Velm = VectorElement('Lagrange', cell, 3)
Qelm = FiniteElement('Lagrange', cell, 2) 
V    = FunctionSpace(sas, Velm)
Q    = FunctionSpace(sas, Qelm) 

# read the sas flow velocity 
vel_file = "results/csf_flow/T1/csf_v.xdmf"
with XDMFFile(vel_file) as file:
    velocity_sas_ = Function(V)
    file.read_checkpoint(velocity_sas_, "velocity")

VP1_elem = VectorElement('Lagrange', cell, 1)
VP1 = FunctionSpace(sas, VP1_elem)
velocity_sas = Function(VP1) 
velocity_sas = interpolate(velocity_sas_ , VP1)

# read the arterial network and mesh 
pvs_ratio_artery  = 2.0 
artery, artery_radii, artery_roots = read_vtk_network("mesh/networks/arteries_smooth.vtk")
artery_radii = as_P0_function(artery_radii)
pvs_radii = Function(artery_radii.function_space())
pvs_radii.vector().set_local(pvs_ratio_artery*artery_radii.vector().get_local())

# define velocity average using Fenics_ii average functionality 
disk_pvs  = Disk(radius = pvs_radii, degree =11)
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

with XDMFFile('results/pvs_flow/uhat_prod.xdmf') as xdmf:
    xdmf.write_checkpoint(pvs_flow_vec, "velocity")

File('results/pvs_flow/uhat_prod.pvd') << pvs_flow_vec
