from dolfin import *
from petsc4py import PETSc
from solver import mark_internal_interface, mark_external_boundary, as_P0_function
from sas_flow import map_on_global
from IPython import embed 
import numpy as np
import xii 

parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['optimize'] = True
parameters["ghost_mode"] = "shared_facet"

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

inlet_id = 2
par_csf_id = 3
par_outer_id = 4
# Parameters
D = Constant(3.8e-10)
t_end = 3600 * 4
dt = 30
beta = Constant(3.8e-7)

mesh = Mesh()

with XDMFFile(f'mesh/T1/volmesh/mesh.xdmf') as f:
    f.read(mesh)
    gdim = mesh.geometric_dimension()
    label = MeshFunction('size_t', mesh, gdim, 0)
    f.read(label, 'label')

with XDMFFile(f'results/csf_flow/sas_flow/csf_v.xdmf') as f:
    V = VectorFunctionSpace(mesh, "CG", 3)
    b = Function(V)
    f.read_checkpoint(b, "velocity")


sm = xii.EmbeddedMesh(label, [CSFID, LVID, V34ID])

with XDMFFile(f'results/csf_flow/cardiac_sas_flow/csf_vis_p.xdmf') as f:
    V = FunctionSpace(sm, "CG", 2)
    p = Function(V)
    f.read_checkpoint(p, "pressure")

    gradp = sqrt(inner(grad(p), grad(p)))
    rho = 993 # kg/m^3
    nu = 7e-7 # m^2/s
    omega = 2*np.pi
    h = 3e-3 / 2
    P = gradp /(rho*omega*nu/h)
    R = P**2
    Sc = nu / D
    alpha = np.sqrt(h**2 * omega / nu)
    beta2 = alpha**2 * Sc

    smDG = FunctionSpace(sm, "CG", 1)
    uR = TrialFunction(smDG)
    vR = TestFunction(smDG)
    aR = (Constant(1e-5)*inner(grad(uR), grad(vR)) + Constant(1)*uR*vR)*dx
    LR = Constant(1)*R*vR*dx
    solR = Function(smDG)
    solR.rename("R","R")
    solve(aR==LR, solR)
    print(assemble(R*dx))
    print(assemble(solR*dx))
    solR = map_on_global(interpolate(solR, FunctionSpace(sm, "DG", 0)), mesh)
    solR.rename("R", "R")
    File("R.pvd") << solR

    D *= (1 + solR)

# Define function spaces
DG = FunctionSpace(mesh, "DG", 1)
v = TestFunction(DG)
u = TrialFunction(DG)

u0 = Function(DG)

# STABILIZATION
h = CellDiameter(mesh)
n = FacetNormal(mesh)
alpha = Constant(1e3)


inlet = CompiledSubDomain("on_boundary && x[2] < zmin + eps",
                           zmin=mesh.coordinates()[:,2].min(), eps=0.4e-3)
bm = MeshFunction("size_t", mesh, 2, 0)
inlet.mark(bm, inlet_id)

for i in [CSFID, LVID, V34ID, CSFNOFLOWID]:
    mark_internal_interface(mesh, label, bm, par_csf_id,
                           doms=[i, PARID])
     
mark_external_boundary(mesh, label, bm, par_outer_id, doms=[PARID])

File("bm.pvd") << bm


ds = Measure("ds", mesh, subdomain_data=bm)
dS = Measure("dS", mesh, subdomain_data=bm)
dx = Measure("dx", mesh, subdomain_data=label)

assert np.isclose(assemble(div(b)*div(b)*dx(CSFID)), 0)
assert np.isclose(assemble(div(b)*div(b)*dx(LVID)), 0)
assert np.isclose(assemble(div(b)*div(b)*dx(V34ID)), 0)

g = Expression(" (t < t1) ? \
              2*c_tot / (t1*t2) * t / A : \
              2*c_tot / (t1*t2) * max(0.0, t2 - t) / A ", 
              t1=3600, t2=7200, c_tot=0.5e-3, t=0, A=assemble(1*ds(inlet_id)), degree=1)

eta = 1.0 
def a(u,v) :
   
    # Bilinear form
    dSi = dS(0)
    a_int = dot(grad(v), D*grad(u))*dx \
        - dot(grad(v),b*u)*dx(CSFID) \
        - dot(grad(v),b*u)*dx(LVID) \
        - dot(grad(v),b*u)*dx(V34ID)

    
    DF = Constant(2)*D('+')*D('-')/(D('+') + D('-'))

    #wavg = lambda f, w: f("+")*w("-")/(w("-") + w("+")) + f("-")*w("+")/(w("-") + w("+"))
    wavg = lambda gr, k: 2*k("+")*k("-") / (k("+") + k("-")) * avg(gr)

    bn = (dot(b, n) + abs(dot(b, n)))/2.0

    a_fac = (alpha/avg(h))*DF*dot(jump(u, n), jump(v, n))*dSi \
            - dot(wavg(grad(u), D), jump(v, n))*dSi \
            - dot(jump(u, n), wavg(grad(v), D))*dSi \
            + beta*jump(u)*jump(v)*dS(par_csf_id)
    

    #a_vel = dot(jump(v), bn('+')*u('+') - bn('-')*u('-') )*dSi
    a_vel = dot(dot(avg(b),n('+'))*avg(u), jump(v))*dSi \
        + (eta/2)*dot(abs(dot(avg(b),n('+')))*jump(u), jump(v))*dSi #\
        #+ dot(v, bn*u)*ds

    a = a_int + a_fac + a_vel

    return a

# Define variational forms

L = (1/dt)*inner(u, v)*dx + a(u,v)
b = (1/dt)*inner(u0, v)*dx + v*g*ds(inlet_id)
 
# Create files for storing results
file = XDMFFile("temp/adv_diff_brain.xdmf")
u = Function(DG)
file.write_checkpoint(u, "c", 0.0)

u.assign(u0)
u.rename("u", "u")

opts = PETSc.Options() 
opts.setValue('ksp_type', 'preonly')    
opts.setValue('pc_type', 'lu')
opts.setValue("pc_factor_mat_solver_type", "mumps")
opts.setValue("mat_mumps_icntl_4", "3")
opts.setValue("mat_mumps_icntl_35", 1)
opts.setValue("mat_mumps_cntl_7",  1e-12)  # BLR eps

A, bb = assemble_system(L, b)
ksp = PETSc.KSP().create()
ksp.setOperators(as_backend_type(A).mat())
ksp.setFromOptions()


# Time-stepping
t = 0.0
i = 0
while t < t_end:
    bb = assemble(b)
    print(t)

    # Compute
    ksp.solve(as_backend_type(bb).vec(), 
              as_backend_type(u.vector()).vec())

    # Move to next time step 
    u0.assign(u)
    t += dt
    g.t = t
    i += 1
    print(assemble(u*dx))
    if i%2==0:
        print(t)
        file.write_checkpoint(u, "c", t, append=True)

file.close()