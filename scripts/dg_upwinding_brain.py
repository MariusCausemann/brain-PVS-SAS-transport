from fenics import *
from ufl import min_value, replace

mesh = Mesh()
with XDMFFile(f'results/csf_flow/sas_flow/csf_v.xdmf') as f:
    f.read(mesh)
    V = VectorFunctionSpace(mesh, "CG", 3)
    u = Function(V)
    f.read_checkpoint(u, "velocity")


u = Constant([0,0, 1])

bcid = 0
inlet_id = 2

inlet = CompiledSubDomain("on_boundary && x[2] < zmin + eps",
                            zmin=mesh.coordinates()[:,2].min(), eps=1e-3)
bm = MeshFunction("size_t", mesh, 2, bcid)
inlet.mark(bm, inlet_id)

V = FunctionSpace(mesh, "DG", 0)

q = project(Constant(0.0), V)
q_init = Function(V).assign(q)

T = 0.1
dt = T/500.0
dtc = Constant(dt)
q_bc = Constant(0.0)
q_in = Constant(1.0)

dq_trial = TrialFunction(V)
phi = TestFunction(V)
a = phi*dq_trial*dx

n = FacetNormal(mesh)
un = 0.5*(dot(u, n) + abs(dot(u, n)))
ds = Measure("ds", mesh, subdomain_data=bm)

L1 = dtc*(q*div(phi*u)*dx
          - conditional(dot(u, n) < 0, phi*dot(u, n)*q_in, 0.0)*ds(inlet_id)
          - conditional(dot(u, n) < 0, phi*dot(u, n)*q_bc, 0.0)*ds(bcid)
          - conditional(dot(u, n) > 0, phi*dot(u, n)*q, 0.0)*ds
          - (phi('+') - phi('-'))*(un('+')*q('+') - un('-')*q('-'))*dS)

q1 = Function(V); q2 = Function(V)
L2 = replace(L1, {q: q1}); L3 = replace(L1, {q: q2})

dq = Function(V)

#params = {'ksp_type': 'preonly', 'pc_type': 'bjacobi', 'sub_pc_type': 'ilu'}
prob1 = LinearVariationalProblem(a, L1, dq)
solv1 = LinearVariationalSolver(prob1)
prob2 = LinearVariationalProblem(a, L2, dq)
solv2 = LinearVariationalSolver(prob2)
prob3 = LinearVariationalProblem(a, L3, dq)
solv3 = LinearVariationalSolver(prob3)

for s in [solv1, solv2, solv3]:
    s.parameters['linear_solver'] = 'mumps'

t = 0.0
step = 0
output_freq = 20
outfile = File("temp/dg_upwinding_brain.pvd")
q.rename("q", "q")
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