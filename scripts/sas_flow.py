from dolfin import *

from solver import * 
import os
import numpy as np
from xii import *

sas = Mesh()
with XDMFFile('mesh/volmesh/mesh.xdmf') as f:
    f.read(sas)
    gdim = sas.geometric_dimension()
    vol_subdomains = MeshFunction('size_t', sas, gdim, 0)
    f.read(vol_subdomains, 'label')
    sas.scale(1e-3)  # scale from mm to m



sas_outer = EmbeddedMesh(vol_subdomains, 1) 



boundary_markers = MeshFunction("size_t", sas, sas.topology().dim() - 1)
boundary_markers.set_all(0)

# Sub domain for efflux route (mark whole boundary of the full domain) 
class Efflux(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

efflux = Efflux()
efflux.mark(boundary_markers,1)

boundary_markers_outer = MeshFunction("size_t", sas_outer, sas_outer.topology().dim() - 1, 0)
DomainBoundary().mark(boundary_markers_outer, 2)
sas_outer.translate_markers(boundary_markers, (1, ), marker_f=boundary_markers_outer)

File("fpp.pvd") << boundary_markers_outer


ds = Measure("ds", domain=sas_outer, subdomain_data=boundary_markers_outer)
# Define function space for velocity and pressure
V = VectorFunctionSpace(sas_outer, 'CG', 2)
Q = FunctionSpace(sas_outer, 'CG', 1)

W = [V, Q]

u , p = map(TrialFunction, W) 
v , q = map(TestFunction,  W)
n = FacetNormal(sas_outer)

mu = Constant(1e-3) # units need to be checked 
R = Constant(-1e-5)
a = block_form(W,2)  
a[0][0] = 2*mu*inner(sym(grad(u)), sym(grad(v)))*dx + R*inner(dot(u,n), dot(v,n))*ds(1)
a[0][1] = -inner(p, div(v))*dx
a[1][0] = -inner(q, div(u))*dx
a[1][1] = Constant(0.0)*inner(p,q)*dx 
L = block_form(W, 1)
# Volumetric
f = Constant((10,1,1)) 
L[0] = inner(f, v)*dx 

V_bcs = [DirichletBC(V, Constant((0.0,0.0,0.0)), ds.subdomain_data(), 2)]
W_bcs = [V_bcs, [] ]
A, b = map(ii_assemble, (a, L))
A, b = apply_bc(A, b, W_bcs)
A, b = map(ii_convert, (A, b))
 
wh = ii_Function(W)


# preconditioner 

# Form for use in constructing preconditioner matrix (H^1 inner product)
P       = block_form(W,2) 
P[0][0] = 2*mu*inner(sym(grad(u)), sym(grad(v)))*dx 
P[1][1] = inner(p,q)*dx 
P, bb  = map(ii_assemble, (P,L)) 
P, bb  = apply_bc(P, bb, W_bcs) 
P ,bb  = map(ii_convert, (P,bb))


krylov_method = "minres"        
solver = KrylovSolver(krylov_method, "amg")
ksp_params = solver.parameters
ksp_params["monitor_convergence"] = True
ksp_params["relative_tolerance"] = 1E-8

solver.set_operators(A, P)
solver.solve(wh.vector(), b) 
# Get sub-functions
u , p = wh 
# Save solutions 
ufile_pvd = File("velocity_sas.pvd")
ufile_pvd << u
pfile_pvd = File("pressure.pvd")
pfile_pvd << p