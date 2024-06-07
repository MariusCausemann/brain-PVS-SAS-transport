import time
import os
import math
from dolfin import *

# Parameters
D = Constant(3.8e-10)
t_end = 1.2
dt = 0.1

# Create mesh and define function space
mesh = UnitSquareMesh(40, 40)
mesh.rotate(10)

# Define function spaces
V = FunctionSpace(mesh, "DG", 1)

b = Constant([0,1])
	
# Define unknown and test function(s)
v = TestFunction(V)
u = TrialFunction(V)

u0 = Function(V)

# STABILIZATION
h = CellDiameter(mesh)
n = FacetNormal(mesh)
alpha = Constant(1)

# ( dot(v, n) + |dot(v, n)| )/2.0
bn = (dot(b, n) + abs(dot(b, n)))/2.0

inlet_id = 2

inlet = CompiledSubDomain("on_boundary && x[0] > 0 && x[1] < 0.1")
bm = MeshFunction("size_t", mesh, 1, 0)
inlet.mark(bm, inlet_id)
ds = Measure("ds", mesh, subdomain_data=bm)
#g = Expression("(x[0] > 0.5 - w/2 && x[0] < 0.5 + w/2) ? -sin((x[0] - 0.5 - w/2) * pi/w) : 0", degree=2, w=0.3)
g = Expression("(x[0] > 0.5 - w/2 && x[0] < 0.5 + w/2) ? 1 : 0", degree=2, w=0.3)

def a(u,v) :
    # Bilinear form
    a_int = dot(grad(v), D*grad(u) - b*u)*dx
    
    upwindgrad = lambda d: conditional(bn("+") > 0, grad(d)("+"), grad(d)("-"))

    a_fac = D*(alpha/avg(h))*dot(jump(u, n), jump(v, n))*dS \
            - D*dot(upwindgrad(u), jump(v, n))*dS \
            - D*dot(jump(u, n), upwindgrad(v))*dS
    
    a_vel = dot(jump(v), bn('+')*u('+') - bn('-')*u('-') )*dS  + dot(v, bn*u)*ds

    a = a_int + a_fac + a_vel
    return a

# Define variational forms
A = (1/dt)*inner(u, v)*dx + a(u,v)
b = (1/dt)*inner(u0,v)*dx + g* v*ds(inlet_id)

# Create files for storing results
file = File("temp/u.pvd")

u = Function(V)
problem = LinearVariationalProblem(A, b, u)
solver  = LinearVariationalSolver(problem)
solver.parameters['linear_solver'] = 'mumps'

u.assign(u0)
u.rename("u", "u")

# Time-stepping
t = 0.0

while t < t_end:

    # Compute
    solver.solve()
    # Save to file
    # Move to next time step
    u0.assign(u)
    t += dt
    file.write(u, t)
    print(u.vector().min())

