from dolfin import *
from petsc4py import PETSc
import ufl 
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['optimize'] = True
parameters["ghost_mode"] = "shared_facet"

# Parameters
D = Constant(3.8e-10)

dt = 100
t_end = 4000 #86400

mesh = Mesh()
with XDMFFile(f'results/csf_flow/sas_flow/csf_v.xdmf') as f:
    f.read(mesh)
    V = VectorFunctionSpace(mesh, "CG", 3)
    b = Function(V)
    f.read_checkpoint(b, "velocity")

#b = Constant([0,0, 1e-5])


# Define function spaces
DG = FunctionSpace(mesh, "DG", 1)
v = TestFunction(DG)
u = TrialFunction(DG)

with XDMFFile(f'temp/u.xdmf') as f:
    #f.read(mesh)
    DG0 = FunctionSpace(mesh, "DG", 0)
    u0 = Function(DG0)
    f.read_checkpoint(u0, "initial")

u0 = interpolate(u0, DG)
#u0 = Function(DG)

# STABILIZATION
h = CellDiameter(mesh)
n = FacetNormal(mesh)
alpha = Constant(200)

theta = Constant(1.0)

# ( dot(v, n) + |dot(v, n)| )/2.0
bn = (dot(b, n) + abs(dot(b, n)))/2.0

bcid = 0
inlet_id = 2
inlet = CompiledSubDomain("on_boundary && x[2] < zmin + eps",
                            zmin=mesh.coordinates()[:,2].min(), eps=1e-3)
bm = MeshFunction("size_t", mesh, 2, bcid)
inlet.mark(bm, inlet_id)
ds = Measure("ds", mesh, subdomain_data=bm)

g = Expression(" (t < t1) ? \
              2*c_tot / (t1*t2) * t / A : \
              2*c_tot / (t1*t2) * max(0.0, t2 - t) / A ", 
              t1=3600, t2=7200, c_tot=0.5e-3, t=0, A=assemble(1*ds(inlet_id)), degree=1)


eta = 1
def a(u,v) :
   
    # Bilinear form
    a_int = dot(grad(v), D*grad(u) - b*u)*dx  - dot(div(b)*u, v)*dx 

    a_fac = D*(alpha/avg(h))*dot(jump(u, n), jump(v, n))*dS \
            - D*dot(avg(grad(u)), jump(v, n))*dS \
            - D*dot(jump(u, n), avg(grad(v)))*dS
        
    #a_vel = dot(jump(v), bn('+')*u('+') - bn('-')*u('-') )*dS  + dot(v, bn*u)*ds
    a_vel = dot(dot(b,n('+'))*avg(u), jump(v))*dS + (eta/2)*dot(abs(dot(b,n('+')))*jump(u), jump(v))*dS + dot(v, bn*u)*ds

    a = a_int + a_fac + a_vel

    return a

# Define variational forms
a0=a(u0,v)
a1=a(u,v)

 
# Create files for storing results
file = File("temp/adv_diff_brain.pvd")
A = (1/dt)*inner(u, v)*dx - (1/dt)*inner(u0,v)*dx + theta*a1 + (1-theta)*a0

F = A  - g* v*ds(inlet_id)

u = Function(DG)

u.assign(u0)
u.rename("u", "u")

opts = PETSc.Options() 
opts.setValue('ksp_type', 'preonly')    
opts.setValue('pc_type', 'lu')
opts.setValue("pc_factor_mat_solver_type", "mumps")
opts.setValue("mat_mumps_icntl_4", "3")

L, f = lhs(F),rhs(F)

A, b = assemble_system(L, f)
ksp = PETSc.KSP().create()
ksp.setOperators(as_backend_type(A).mat())
ksp.setFromOptions()

# Time-stepping
t = 0.0

file << u

i = 0
while t < t_end:
    b = assemble(f)
    # Compute
    ksp.solve(as_backend_type(b).vec(), 
              as_backend_type(u.vector()).vec())

    # Move to next time step
    u0.assign(u)
    t += dt
    g.t = t
    i += 1
    if i%10==0:
        # Save to file
        print(u.vector().min())
        file.write(u, t)


#with XDMFFile(f'temp/u.xdmf') as xdmf:
#    xdmf.write_checkpoint(u, "initial")