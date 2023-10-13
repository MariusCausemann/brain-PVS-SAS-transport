from petsc4py import PETSc
from dolfin import *
import xii 
from solver import * 
import os

if __name__ == '__main__':
   
    vein, vein_radii, vein_roots = read_vtk_network("../mesh/networks/venes_smooth.vtk")
    vein_radii = as_P0_function(vein_radii)

    artery, artery_radii, artery_roots = read_vtk_network("../mesh/networks/arteries_smooth.vtk")
    artery_radii = as_P0_function(artery_radii)

    sas = Mesh()
    with XDMFFile('../mesh/volmesh/mesh.xdmf') as f:
        f.read(sas)
        sas_subdomains = MeshFunction('size_t', sas, 3,0)
        f.read(sas_subdomains, 'label')

    Ds = pcws_constant(sas_subdomains, {1: Constant(20),  # sas
                                        2: Constant(30)})

    V = FunctionSpace(sas, 'CG', 1)
    Qa = FunctionSpace(artery, 'CG',1)
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
    dt = 0.001 
    num_timesteps = 300
    T  = num_timesteps*dt  
    a = xii.block_form(W, 2)
    # Ds = Constant(2.0) 
    Da = Constant(1.0) 
    Dv = Constant(1.0)
    velocity_a = Constant(1.0)
    velocity_v = Constant(1.0)
    fa = Constant(10.0) 
    fv = Constant(10.0)

    hF, nF, gamma_Nitsche = CellDiameter(sas), FacetNormal(sas), Constant(100)
    gSAS = Constant(0.0)
    xi = Constant(17.0)
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

    
    V_bcs  =  [DirichletBC(V,Constant(0.0), DomainBoundary())]
    Qa_bcs = []
    Qv_bcs = []
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
    opts.setValue('ksp_rtol', 1E-40)
    opts.setValue('ksp_atol', 1E-12)   # |AX-b| < 1E-12
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
    vtkfile_4 = File(results_dir + 'uh_real_sas.pvd') 
    vtkfile_5 = File(results_dir + 'uh_hemi1.pvd') 
    #vtkfile_6 = File('uh_hemi2.pvd') 

    u_i.rename("SAS", "time")
    pa_i.rename("Arteries", "time")
    pv_i.rename("Veins", "time")

    vtkfile_1 << (u_i, float(0.0))
    vtkfile_2 << (pa_i, float(0.0))
    vtkfile_3 << (pv_i, float(0.0))

    sas_mesh = SubMesh(sas, sas_subdomains, 1)
    V_sas = FunctionSpace(sas_mesh, V.ufl_element())
    hemi1_mesh = SubMesh(sas, sas_subdomains,2) 
    V_hemi1    = FunctionSpace(hemi1_mesh, V.ufl_element())
    #hemi2_mesh = SubMesh(sas, sas_subdomains,3) 
    #V_hemi2    = FunctionSpace(hemi2_mesh, V.ufl_element())
    uh_sas = Function(V_sas)
    uh_sas.rename('SAS', 'time')
    uh_hemi1 = Function(V_hemi1)
    uh_hemi1.rename('Hemispheres', 'time')

    #uh_hemi2 = Function(V_hemi2)
   # uh_hemi2.rename('Hemispheres', 'time')

    uh_hemi1.assign(interpolate(u_i, V_hemi1)) 
    #uh_hemi2.assign(interpolate(u_i, V_hemi2)) 
    vtkfile_5 << (uh_hemi1, 0.0 )
   # vtkfile_6 << (uh_hemi2, 0.0)
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

        uh_sas.assign(interpolate(wh[0], V_sas))
        uh_hemi1.assign(interpolate(wh[0], V_hemi1))
        #uh_hemi2.assign(interpolate(wh[0], V_hemi2))

        t += dt 
        wh[0].rename("SAS", "time")
        wh[1].rename("Arteries", "time")
        wh[2].rename("Veins", "time")
        vtkfile_1 << (wh[0],float(t)) 
        vtkfile_2 << (wh[1],float(t)) 
        vtkfile_3 << (wh[2],float(t)) 
        vtkfile_4 << (uh_sas,float(t)) 
        vtkfile_5 << (uh_hemi1,float(t)) 
        #vtkfile_6 << (uh_hemi2,float(t))    