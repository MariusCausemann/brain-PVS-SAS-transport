from dolfin import *

from solver import * 
import os
import numpy as np
from xii import *
from petsc4py import PETSc
import sympy as sp



# get mesh 
sas = Mesh()
with XDMFFile('mesh/volmesh/mesh.xdmf') as f:
    f.read(sas)
    gdim = sas.geometric_dimension()
    vol_subdomains = MeshFunction('size_t', sas, gdim, 0)
    f.read(vol_subdomains, 'label')
    sas.scale(1e-3)  # scale from mm to m

# pick the sas mesh 
sas_outer = EmbeddedMesh(vol_subdomains, 1) 

# create boundary markers 
boundary_markers = MeshFunction("size_t", sas, sas.topology().dim() - 1)
boundary_markers.set_all(0)
# Sub domain for efflux route (mark whole boundary of the full domain) 
class Efflux(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary
efflux = Efflux()
efflux.mark(boundary_markers,1)

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
Welm = MixedElement([Velm, Qelm]) 
W = FunctionSpace(sas_outer, Welm)

u , p = TrialFunctions(W) 
v , q = TestFunctions(W)
n = FacetNormal(sas_outer)

mu = Constant(0.7e-3) # units need to be checked 
R = Constant(-1e-5)
f = Constant((0,0,0)) # right hand side, ofcourse needs to modified to obtain the right thing 
cp1_midpoint = [0.128, 0.229, 0.192] # found in paraview
cp2_midpoint = [0.2, 0.229, 0.192] # found in paraview

A = 0.0007475317873778963 # normalize to reach 0.5 L / day
g1 =  Expression("A*exp( - ((x[0] - m0)*(x[0] - m0) + (x[1] - m1)*(x[1] - m1) + (x[2] - m2)*(x[2] - m2)) / (sigma*sigma))",
                    m0=cp1_midpoint[0], m1=cp1_midpoint[1],
                    m2=cp1_midpoint[2], sigma=0.01, A=A, degree=2)
g2 =  Expression("A*exp( - ((x[0] - m0)*(x[0] - m0) + (x[1] - m1)*(x[1] - m1) + (x[2] - m2)*(x[2] - m2)) / (sigma*sigma))",
                    m0=cp2_midpoint[0], m1=cp2_midpoint[1],
                    m2=cp2_midpoint[2], sigma=0.01, A=A, degree=2)
g = g1 + g2
gtot = assemble(g*dx(domain=sas_outer))
g_L_per_day = gtot * 1e3 *(60*60*24)
print(f"production rate: {g_L_per_day} L/day")
assert np.isclose(g_L_per_day, 0.5)
gfunc = project(g, Q)
File("g.pvd") << gfunc
a = (inner(2*mu*sym(grad(u)), sym(grad(v)))*dx - inner(p, div(v))*dx
         -inner(q, div(u))*dx  + inner(R*dot(u,n), dot(v,n))*ds(1)) 
bcs = [DirichletBC(W.sub(0), Constant((0.0,0.0,0.0)), ds.subdomain_data(), 2)] 

L = inner(f, v)*dx + g*q*dx

A, b = assemble_system(a, L, bcs)

# Preconditioner
a_prec = (inner(2*mu*grad(u), grad(v))*dx + inner(R*dot(u,n), dot(v,n))*ds(1)
            + (1/mu)*inner(p, q)*dx) 

B, _ = assemble_system(a_prec, L, bcs)    


print(W.sub(0).collapse().dim(), W.sub(1).collapse().dim())

wh = Function(W)

solver = PETScKrylovSolver()
solver.set_operators(A, B)

ksp = solver.ksp()

opts = PETSc.Options()
opts.setValue('ksp_type', 'fgmres')
opts.setValue('ksp_atol', 1E-10)       
opts.setValue('ksp_rtol', 1E-40)            
opts.setValue('ksp_view_pre', None)
opts.setValue('ksp_norm_type', 'unpreconditioned')
opts.setValue('ksp_monitor_true_residual', None)                
opts.setValue('ksp_converged_reason', None)
# Specify that the preconditioner is block diagonal and customize
# inverses of individual blocks
opts.setValue('fieldsplit_0_ksp_type', 'preonly')
#opts.setValue('fieldsplit_0_pc_type', 'lu')        # Exact inverse for velocity
opts.setValue('fieldsplit_0_pc_type', 'hypre')        # Exact inverse for velocity

#   -pc_hypre_boomeramg_cycle_type <V> (choose one of) V W (None)
#   -pc_hypre_boomeramg_max_levels <25>: Number of levels (of grids) allowed (None)
#opts.setValue('fieldsplit_0_pc_hypre_boomeramg_max_iter', 2)  # : Maximum iterations used PER hypre call (None)
#   -pc_hypre_boomeramg_tol <0.>: Convergence tolerance PER hypre call (0.0 = use a fixed number of iterations) (None)
#   -pc_hypre_boomeramg_truncfactor <0.>: Truncation factor for interpolation (0=no truncation) (None)
#   -pc_hypre_boomeramg_P_max <0>: Max elements per row for interpolation operator (0=unlimited) (None)
#   -pc_hypre_boomeramg_agg_nl <0>: Number of levels of aggressive coarsening (None)
#   -pc_hypre_boomeramg_agg_num_paths <1>: Number of paths for aggressive coarsening (None)
#opts.setValue('fieldsplit_0_pc_hypre_boomeramg_strong_threshold', 0.7)  # : Threshold for being strongly connected (None)
#   -pc_hypre_boomeramg_max_row_sum <0.9>: Maximum row sum (None)
#   -pc_hypre_boomeramg_grid_sweeps_all <1>: Number of sweeps for the up and down grid levels (None)
#opts.setValue('fieldsplit_0_pc_hypre_boomeramg_nodal_coarsen', 1) #  Use a nodal based coarsening 1-6 (HYPRE_BoomerAMGSetNodal)
#   -pc_hypre_boomeramg_vec_interp_variant <0>: Variant of algorithm 1-3 (HYPRE_BoomerAMGSetInterpVecVariant)
#   -pc_hypre_boomeramg_grid_sweeps_down <1>: Number of sweeps for the down cycles (None)\
#opts.setValue('fieldsplit_0_pc_hypre_boomeramg_grid_sweeps_down', 3)
#opts.setValue('fieldsplit_0_pc_hypre_boomeramg_grid_sweeps_up', 3) # Number of sweeps for the up cycles (None)
#   -pc_hypre_boomeramg_grid_sweeps_coarse <1>: Number of sweeps for the coarse level (None)
#   -pc_hypre_boomeramg_smooth_type <Schwarz-smoothers> (choose one of) Schwarz-smoothers Pilut ParaSails Euclid (None)
#   -pc_hypre_boomeramg_smooth_num_levels <25>: Number of levels on which more complex smoothers are used (None)
#   -pc_hypre_boomeramg_eu_level <0>: Number of levels for ILU(k) in Euclid smoother (None)
#   -pc_hypre_boomeramg_eu_droptolerance <0.>: Drop tolerance for ILU(k) in Euclid smoother (None)
#   -pc_hypre_boomeramg_eu_bj: <FALSE> Use Block Jacobi for ILU in Euclid smoother? (None)
#   -pc_hypre_boomeramg_relax_type_all <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
#   -pc_hypre_boomeramg_relax_type_down <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
#   -pc_hypre_boomeramg_relax_type_up <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
#   -pc_hypre_boomeramg_relax_type_coarse <Gaussian-elimination> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
#   -pc_hypre_boomeramg_relax_weight_all <1.>: Relaxation weight for all levels (0 = hypre estimates, -k = determined with k CG steps) (None)
#   -pc_hypre_boomeramg_relax_weight_level <1.>: Set the relaxation weight for a particular level (weight,level) (None)
#   -pc_hypre_boomeramg_outer_relax_weight_all <1.>: Outer relaxation weight for all levels (-k = determined with k CG steps) (None)
#   -pc_hypre_boomeramg_outer_relax_weight_level <1.>: Set the outer relaxation weight for a particular level (weight,level) (None)
#   -pc_hypre_boomeramg_no_CF: <FALSE> Do not use CF-relaxation (None)
#   -pc_hypre_boomeramg_measure_type <local> (choose one of) local global (None)
#   -pc_hypre_boomeramg_coarsen_type <Falgout> (choose one of) CLJP Ruge-Stueben  modifiedRuge-Stueben   Falgout  PMIS  HMIS (None)
#   -pc_hypre_boomeramg_interp_type <classical> (choose one of) classical   direct multipass multipass-wts ext+i ext+i-cc standard standard-wts   FF FF1 (None)
#   -pc_hypre_boomeramg_print_statistics: Print statistics (None)
#   -pc_hypre_boomeramg_print_statistics <3>: Print statistics (None)
#   -pc_hypre_boomeramg_print_debug: Print debug information (None)
#   -pc_hypre_boomeramg_nodal_relaxation: <FALSE> Nodal relaxation via Schwarz (None)

opts.setValue('fieldsplit_1_ksp_type', 'preonly')
opts.setValue('fieldsplit_1_pc_type', 'lu')     # AMG for pressure

#opts.setValue('fieldsplit_2_ksp_type', 'preonly')
#opts.setValue('fieldsplit_2_pc_type', 'lu')                
pc = ksp.getPC()
pc.setType(PETSc.PC.Type.FIELDSPLIT)
is_V = PETSc.IS().createGeneral(W.sub(0).dofmap().dofs())
is_Q = PETSc.IS().createGeneral(W.sub(1).dofmap().dofs())

pc.setFieldSplitIS(('0', is_V), ('1', is_Q))
pc.setFieldSplitType(PETSc.PC.CompositeType.MULTIPLICATIVE) 

ksp.setUp()

subksps = pc.getFieldSplitSubKSP()

A0, B0 = subksps[0].getOperators()

gdim = W.mesh().geometry().dim()
A0.setBlockSize(gdim)
B0.setBlockSize(gdim)

pc.setFromOptions()
ksp.setFromOptions()
niters = solver.solve(wh.vector(), b)
rnorm = ksp.getResidualNorm()
uh, ph = wh.split(deepcopy=True)[:2]
# write solutions 
ufile_pvd = File("velocity_sas.pvd")
ufile_pvd << uh
pfile_pvd = File("pressure.pvd")
pfile_pvd << ph


with XDMFFile('results/sas_flow.xdmf') as xdmf:
    xdmf.write_checkpoint(uh, "velocity")
with XDMFFile('results/sas_flow_vis.xdmf') as xdmf:
    xdmf.write(uh)


with XDMFFile('results/sas_p.xdmf') as xdmf:
    xdmf.write_checkpoint(ph, "pressure")
with XDMFFile('results/sas_p_vis.xdmf') as xdmf:
    xdmf.write(ph)
