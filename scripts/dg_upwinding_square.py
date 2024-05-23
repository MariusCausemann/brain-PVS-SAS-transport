from fenics import *
from ufl import min_value, replace

mesh = UnitSquareMesh(40, 40)

V = FunctionSpace(mesh, "DG", 0)
W = VectorFunctionSpace(mesh, "CG", 1)

x, y = SpatialCoordinate(mesh)

velocity = as_vector((0.5 - y, x - 0.5))
u = project(velocity, W)

cone_r0 = 0.15; cone_x0 = 0.5; cone_y0 = 0.25
cyl_r0 = 0.15; cyl_x0 = 0.5; cyl_y0 = 0.75

cone = 1.0 - min_value(sqrt(pow(x-cone_x0, 2) + pow(y-cone_y0, 2))/cyl_r0, 1.0)

q = project(1.0 + cone, V)
q_init = Function(V).assign(q)

T = 2*pi
dt = T/600.0
dtc = Constant(dt)
q_in = Constant(2.0)

dq_trial = TrialFunction(V)
phi = TestFunction(V)
a = phi*dq_trial*dx

n = FacetNormal(mesh)
un = 0.5*(dot(u, n) + abs(dot(u, n)))

L1 = dtc*(q*div(phi*u)*dx
          - conditional(dot(u, n) < 0, phi*dot(u, n)*q_in, 0.0)*ds
          - conditional(dot(u, n) > 0, phi*dot(u, n)*q, 0.0)*ds
          - (phi('+') - phi('-'))*(un('+')*q('+') - un('-')*q('-'))*dS)

q1 = Function(V); q2 = Function(V)
L2 = replace(L1, {q: q1}); L3 = replace(L1, {q: q2})

dq = Function(V)

params = {'ksp_type': 'preonly', 'pc_type': 'bjacobi', 'sub_pc_type': 'ilu'}
prob1 = LinearVariationalProblem(a, L1, dq)
solv1 = LinearVariationalSolver(prob1)
prob2 = LinearVariationalProblem(a, L2, dq)
solv2 = LinearVariationalSolver(prob2)
prob3 = LinearVariationalProblem(a, L3, dq)
solv3 = LinearVariationalSolver(prob3)

for s in [solv1, solv2, solv3]:
    s.parameters['linear_solver'] = 'lu'

t = 0.0
step = 0
output_freq = 20
outfile = File("temp/dg_upwinding_square.pvd")
while t < T - 0.5*dt:
    solv1.solve()
    q1.assign(q + dq)

    solv2.solve()
    q2.assign(0.75*q + 0.25*(q1 + dq))

    solv3.solve()
    q.assign((1.0/3.0)*q + (2.0/3.0)*(q2 + dq))

    step += 1
    t += dt

    if step % output_freq == 0:
        print("t=", t)
    outfile.write(q, t)
