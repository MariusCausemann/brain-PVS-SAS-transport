from petsc4py import PETSc
from dolfin import *
import xii 
from solver import * 
import os
import numpy as np

if __name__ == '__main__':
   
    vein, vein_radii, vein_roots = read_vtk_network("../mesh/networks/venes_smooth.vtk")
    vein_radii = as_P0_function(vein_radii)

    artery, artery_radii, artery_roots = read_vtk_network("../mesh/networks/arteries_smooth.vtk")
    artery_radii = as_P0_function(artery_radii)

    sas = Mesh()
    with XDMFFile('../mesh/volmesh/mesh.xdmf') as f:
        f.read(sas)
        sas_subdomains = MeshFunction('size_t', sas, 3, 0)
        f.read(sas_subdomains, 'label')

    assert np.allclose(np.unique(sas_subdomains.array()), [1,2])
    [m.scale(1e-3) for m in [sas, vein, artery]] # scale from mm to m

    dt = 60
    T = 60*60*1 # 4h
    num_timesteps = int(T / dt)

    m = (0.17, 0.21, 0.1)
    inlet_radius = 0.040

    inlet = CompiledSubDomain("on_boundary && " + 
                            "+ (x[0] - m0)*(x[0] - m0) " + 
                            "+ (x[1] - m1)*(x[1] - m1) " + 
                            "+ (x[2] - m2)*(x[2] - m2) < r*r",
                            m0=m[0], m1=m[1], m2=m[2], r=inlet_radius)
    
    bm = MeshFunction("size_t", sas, 2, 0)
    inlet.mark(bm, 1)
    
    boundary_concentration = Constant(1)

    mm2m = 1e-3
    ecs_vol = 0.2
    Ds = pcws_constant(sas_subdomains, {1: Constant(3.8*1e-4 * mm2m),  # csf
                                        2: Constant(1.2*1e-4 * mm2m*ecs_vol) # parenchyma
                                        })
        
    Da = Constant(3.8*1e-4 * mm2m) 
    Dv = Constant(3.8*1e-4 * mm2m)
    xi = Constant(1e-3)

    velocity_a = Constant(0.0)
    velocity_v = Constant(0.0)

    fa = Constant(0.0) 
    fv = Constant(0.0)
    gSAS = Constant(0.0)

    V = FunctionSpace(sas, 'CG', 1)
    Qa = FunctionSpace(artery, 'CG', 1)
    Qv = FunctionSpace(vein, 'CG', 1)
    W = [V, Qa, Qv]

    u, pa, pv = map(TrialFunction, W)
    v, qa, qv = map(TestFunction, W)
     
    # initial conditions 
    u_o  = Constant(0.0)
    pa_o = Constant(0.0) 
    pv_o = Constant(0.0) 

    u_i  = interpolate(u_o, V) 
    pa_i = interpolate(pa_o, Qa) 
    pv_i = interpolate(pv_o, Qv)
    # Things for restriction
    dx_a = Measure('dx', domain=artery)
    artery_shape = xii.Circle(radius=artery_radii, degree=20)
    ua, va = (xii.Average(x, artery, artery_shape) for x in (u, v))

    dx_v = Measure('dx', domain=vein)
    vein_shape = xii.Circle(radius=vein_radii, degree=20)
    uv, vv = (xii.Average(x, vein, vein_shape) for x in (u, v))

    a = xii.block_form(W, 2)


    hF, nF, gamma_Nitsche = CellDiameter(sas), FacetNormal(sas), Constant(100)

    a[0][0] = (1/dt)*inner(u,v)*dx + Ds*inner(grad(u), grad(v))*dx + xi*inner(ua, va)*dx_a + xi*inner(uv, vv)*dx_v
    a[0][1] = -xi*inner(pa, va)*dx_a
    a[0][2] = -xi*inner(pv, vv)*dx_v

    a[1][0] = -xi*inner(qa, ua)*dx_a
    a[1][1] = (1/dt)*inner(pa,qa)*dx + Da*inner(grad(pa), grad(qa))*dx - inner(pa, dot(velocity_a,grad(qa)[0]))*dx + xi*inner(pa, qa)*dx

    a[2][0] = -xi*inner(qv, uv)*dx_v
    a[2][2] = (1/dt)*inner(pv,qv)*dx + Dv*inner(grad(pv), grad(qv))*dx - inner(pv, dot(velocity_v,grad(qv)[0]))*dx + xi*inner(pv, qv)*dx 

    L = xii.block_form(W, 1)
    L[0]  = (1/dt)*inner(u_i,v)*dx 
    
    L[1]  = (1/dt)*inner(pa_i, qa)*dx + inner(fa,qa)*dx
    L[2]  = (1/dt)*inner(pv_i, qv)*dx + inner(fv,qv)*dx


    AA, bb = map(xii.ii_assemble, (a, L))

    
    V_bcs  =  [DirichletBC(V, boundary_concentration, inlet)]
    Qa_bcs = [DirichletBC(Qa, boundary_concentration, inlet)]
    Qv_bcs = [DirichletBC(Qv, boundary_concentration, inlet)]
    W_bcs = [V_bcs, Qa_bcs, Qv_bcs]

    AA, _, bc_apply_b = xii.apply_bc(AA, bb, bcs=W_bcs, return_apply_b=True)

    A_ = ksp_mat(xii.ii_convert(AA))

    opts = PETSc.Options()
    opts.setValue('ksp_type', 'cg')    
    #opts.setValue('ksp_view', None)
    #opts.setValue('ksp_view_eigenvalues', None)
    #opts.setValue('ksp_converged_reason', None)
    #opts.setValue('ksp_norm_type', 'unpreconditioned')
    #opts.setValue('ksp_monitor_true_residual', None)
    #opts.setValue('ksp_rtol', 1E-40)
    opts.setValue('ksp_atol', 1E-12)   # |AX-b| < 1E-
    opts.setValue('pc_type', 'hypre')
    opts.setValue('ksp_initial_guess_nonzero', 1)

    ksp = PETSc.KSP().create()
    ksp.setOperators(A_, A_)
    ksp.setFromOptions()
    print('Start solve')
    t = 0.0 

    results_dir = "../results/"
    os.makedirs(results_dir, exist_ok=True)
    vtkfile_1 = File(results_dir + 'uh_sas.pvd') 
    vtkfile_2 = File(results_dir + 'uh_artery.pvd')
    vtkfile_3 = File(results_dir + 'uh_vein.pvd') 

    u_i.rename("c", "time")
    pa_i.rename("c", "time")
    pv_i.rename("c", "time")

    vtkfile_1 << (u_i, float(0.0))
    vtkfile_2 << (pa_i, float(0.0))
    vtkfile_3 << (pv_i, float(0.0))

    wh = xii.ii_Function(W)
    x_ = A_.createVecLeft()
    while t < T: 
        print("time", t)
        bb = xii.ii_assemble(L)
        b = bc_apply_b(bb)

        b_ = ksp_vec(xii.ii_convert(b))
        ksp.solve(b_, x_)
        # NOTE: solve(b_, ksp_vec(wh.vector())) segfault most likely because
        # of type incompatibility seq is expected and we have nest
        wh.vector()[:] = PETScVector(x_)
        u_i.assign(wh[0]) 
        pa_i.assign(wh[1]) 
        pv_i.assign(wh[2])

        t += dt 
        wh[0].rename("c", "time")
        wh[1].rename("c", "time")
        wh[2].rename("c", "time")
        vtkfile_1 << (wh[0],float(t)) 
        vtkfile_2 << (wh[1],float(t)) 
        vtkfile_3 << (wh[2],float(t))