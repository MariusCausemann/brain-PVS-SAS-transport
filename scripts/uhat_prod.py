from dolfin import *

from solver import * 
import os
import numpy as np
from xii import *
from petsc4py import PETSc
import sympy as sp


# get mesh 
sas = Mesh()
with XDMFFile('mesh/mid_mesh/volmesh/mesh.xdmf') as f:
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
vel_file = "results/csf_flow/mid_mesh/csf_v.xdmf"
with XDMFFile(vel_file) as file:
    velocity_sas_ = Function(V)
    file.read_checkpoint(velocity_sas_, "velocity")

VP1_elem = VectorElement('Lagrange', cell, 1)
VP1 = FunctionSpace(sas, VP1_elem)
velocity_sas = Function(VP1) 
velocity_sas = interpolate(velocity_sas_ , VP1)

File('results/velocity_test_sas.pvd') << velocity_sas
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

with XDMFFile('results/pvs_flow/pvs_flow.xdmf') as xdmf:
    xdmf.write_checkpoint(pvs_flow, "velocity")


File('results/pvs_flow/pvs_flow.pvd') << pvs_flow
