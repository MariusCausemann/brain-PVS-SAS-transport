from dolfin import *
from petsc4py import PETSc
from solver import mark_internal_interface, mark_external_boundary

parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['optimize'] = True
parameters["ghost_mode"] = "shared_facet"

# Parameters
D = Constant(3.8e-10)
t_end = 7200
dt = 30
beta = Constant(3.8e-7)

mesh = Mesh()

from IPython import embed 
# embed()


with XDMFFile(f'results/csf_flow/sas_flow/csf_v.xdmf') as f:
    f.read(mesh)
    V = VectorFunctionSpace(mesh, "CR", 1)
    b = Function(V)
    f.read_checkpoint(b, "velocity")
 
with XDMFFile(f'results/csf_flow/cardiac_sas_flow/csf_v.xdmf') as f:
    cb = Function(V)
    f.read_checkpoint(cb, "velocity")


umag = project(sqrt(inner(cb, cb)), FunctionSpace(mesh, "DG", 1),
                solver_type="cg", preconditioner_type="jacobi")

dispersion_fac = 1e3 / umag.vector().max()
D *= (1 + umag*dispersion_fac)

with XDMFFile(f'mesh/T1/volmesh/mesh.xdmf') as f:
    gdim = mesh.geometric_dimension()
    label = MeshFunction('size_t', mesh, gdim, 1)
    f.read(label, 'label')

with XDMFFile(f'mesh_labeled.xdmf') as f: 
    bm = MeshFunction('size_t', mesh, gdim - 1,0) 
    f.read(bm)

File("bm.pvd") << bm
# Define function spaces
DG = FunctionSpace(mesh, "DG", 1)
v = TestFunction(DG)
u = TrialFunction(DG)

u0 = Function(DG)

# STABILIZATION
h = CellDiameter(mesh)
n = FacetNormal(mesh)
alpha = Constant(100)

# ( dot(v, n) + |dot(v, n)| )/2.0
bn = (dot(b, n) + abs(dot(b, n)))/2.0

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

inlet_id = 2
par_csf_id = 3
par_outer_id = 4

#inlet = CompiledSubDomain("on_boundary && x[2] < zmin + eps",
#                           zmin=mesh.coordinates()[:,2].min(), eps=0.4e-3)
#bm = MeshFunction("size_t", mesh, 2, 0)
#inlet.mark(bm, inlet_id)

#for i in [CSFID, LVID, V34ID, CSFNOFLOWID]:
#    mark_internal_interface(mesh, label, bm, par_csf_id,
#                           doms=[i, PARID])
     
#mark_external_boundary(mesh, label, bm, par_outer_id, doms=[PARID])

File("bm.pvd") << bm


ds = Measure("ds", mesh, subdomain_data=bm)
dS = Measure("dS", mesh, subdomain_data=bm)
dx = Measure("dx", mesh, subdomain_data=label)


g = Expression(" (t < t1) ? \
              2*c_tot / (t1*t2) * t / A : \
              2*c_tot / (t1*t2) * max(0.0, t2 - t) / A ", 
              t1=3600, t2=7200, c_tot=0.5e-3, t=0, A=assemble(1*ds(inlet_id)), degree=1)

#g = Expression('exp(- (pow(x[0] - c0, 2) + pow(x[1] - c1, 2))/(2*sig))',
#                      degree=4, sig=0.000002, c0=coords[0], c1=coords[1])
gf = interpolate(g, DG)
gf.rename("g", "g")

eta = 1.0 
def a(u,v) :
   
    # Bilinear form
    dSi = dS(0)
    a_int = dot(grad(v), D*grad(u) - b*u)*dx
    
    DF = Constant(2)*D('+')*D('-')/(D('+') + D('-'))

    wavg = lambda f, w: f("+")*w("-")/(w("-") + w("+")) + f("-")*w("+")/(w("-") + w("+"))

    a_fac = (alpha/avg(h))*dot(jump(u, n), jump(v, n))*dSi \
            - dot(avg(D*grad(u)), jump(v, n))*dSi \
            - dot(jump(u, n), avg(D*grad(v)))*dSi \
            + beta*jump(u)*jump(v)*dS(par_csf_id)
    
    a_vel = dot(jump(v), bn('+')*u('+') - bn('-')*u('-') )*dSi +  dot(v, bn*u)*ds
    #a_vel = dot(dot(b('+'),n('+'))*avg(u), jump(v))*dS + (eta/2)*dot(abs(dot(b('+'),n('+')))*jump(u), jump(v))*dS + dot(v, bn*u)*ds

    a = a_int + a_fac + a_vel

    return a

# Define variational forms

L = (1/dt)*inner(u, v)*dx + a(u,v)
b = (1/dt)*inner(u0, v)*dx + v*g*ds(inlet_id)
 
# Create files for storing results
file = XDMFFile("temp/adv_diff_brain.xdmf")
u = Function(DG)
file.write_checkpoint(u, "velocity", 0.0)

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
    if i%5==0:
        file.write_checkpoint(u, "velocity", t, append=True)

file.close()