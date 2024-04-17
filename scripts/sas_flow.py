from dolfin import *

from solver import * 
import os
import numpy as np
from xii import *
from petsc4py import PETSc
import sympy as sp

PETScOptions.set("mat_mumps_icntl_4", 3)  # mumps verbosity
PETScOptions.set("mat_mumps_icntl_24", 1)  # null pivot detection
#PETScOptions.set("mat_mumps_icntl_11", 1)  # error analysis

# get mesh 
sas = Mesh()
with XDMFFile('mesh/mid_mesh/volmesh/mesh.xdmf') as f:
    f.read(sas)
    gdim = sas.geometric_dimension()
    sas_components = MeshFunction('size_t', sas, gdim, 0)
    f.read(sas_components, 'sas_components')
    sas.scale(1e-3)  # scale from mm to m

# pick the sas mesh 
sas_outer = EmbeddedMesh(sas_components, 1) 

# create boundary markers 
boundary_markers = MeshFunction("size_t", sas, sas.topology().dim() - 1)
boundary_markers.set_all(0)
# Sub domain for efflux route (mark whole boundary of the full domain) 
class Efflux(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[2] > 0.2
efflux = Efflux()
efflux.mark(boundary_markers, 1)

# translate markers to the sas outer mesh 
boundary_markers_outer = MeshFunction("size_t", sas_outer, sas_outer.topology().dim() - 1, 0)
DomainBoundary().mark(boundary_markers_outer, 2)
sas_outer.translate_markers(boundary_markers, (1, ), marker_f=boundary_markers_outer)

File("fpp.pvd") << boundary_markers_outer

ds = Measure("ds", domain=sas_outer, subdomain_data=boundary_markers_outer)


# Define function spaces for velocity and pressure
cell = sas_outer.ufl_cell()  
Velm = VectorElement('Lagrange', cell, 2)
Qelm = FiniteElement('Lagrange', cell, 1) 
Q = FunctionSpace(sas_outer, Qelm)
W = FunctionSpace(sas_outer, Velm * Qelm)

u , p = TrialFunctions(W) 
v , q = TestFunctions(W)
n = FacetNormal(sas_outer)

mu = Constant(0.7e-3) # units need to be checked 
R = Constant(-1e-5)
cp1_midpoint = [0.128, 0.229, 0.192] # found in paraview
cp2_midpoint = [0.2, 0.229, 0.192] # found in paraview

g1 =  Expression("exp( - ((x[0] - m0)*(x[0] - m0) + (x[1] - m1)*(x[1] - m1) + (x[2] - m2)*(x[2] - m2)) / (sigma*sigma))",
                    m0=cp1_midpoint[0], m1=cp1_midpoint[1],
                    m2=cp1_midpoint[2], sigma=0.01, degree=3)
g2 =  Expression("exp( - ((x[0] - m0)*(x[0] - m0) + (x[1] - m1)*(x[1] - m1) + (x[2] - m2)*(x[2] - m2)) / (sigma*sigma))",
                    m0=cp2_midpoint[0], m1=cp2_midpoint[1],
                    m2=cp2_midpoint[2], sigma=0.01, degree=3)


gtot = assemble((g1 + g2 )*dx(domain=sas_outer))

total_production = 0.5e-3/(60*60*24) # 0.5 L / day
g = (g1+ g2)* Constant(total_production / gtot)

f = Constant([0]*gdim)
g_L_per_day = 1e3 *(60*60*24) * assemble(g*dx(domain=sas_outer))
print(f"production rate: {g_L_per_day} L/day")
assert np.isclose(g_L_per_day, 0.5)
gfunc = project(g, Q)
File("g.pvd") << gfunc
a = (inner(2*mu*sym(grad(u)), sym(grad(v)))*dx - inner(p, div(v))*dx
         -inner(q, div(u))*dx + inner(R*dot(u,n), dot(v,n))*ds(1)) 
L = inner(f, v)*dx - g*q*dx
bcs = [DirichletBC(W.sub(0), Constant((0, 0, 0)), ds.subdomain_data(), 2)] 

wh = Function(W)

solve(a==L, wh, bcs=bcs, solver_parameters={"linear_solver":"mumps"})

uh, ph = wh.split(deepcopy=True)[:2]

with XDMFFile('results/sas_flow.xdmf') as xdmf:
    xdmf.write_checkpoint(uh, "velocity")

with XDMFFile('results/sas_p.xdmf') as xdmf:
    xdmf.write_checkpoint(ph, "pressure")
