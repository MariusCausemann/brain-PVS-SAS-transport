from petsc4py import PETSc
from dolfin import *
import xii 
from solver import * 
import os
import numpy as np
import typer
from pathlib import Path
from plotting_utils import read_config

def run_simulation(configfile: str):

    config = read_config(configfile)
    modelname = Path(configfile).stem

    dt = config["dt"]
    T = config["T"]

    inlet_midpoint = config["inlet_midpoint"]
    inlet_radius = config["inlet_radius"]

    ecs_share = config["ecs_share"]
    sas_diffusion = config["sas_diffusion"]
    arterial_pvs_diffusion = config["arterial_pvs_diffusion"]
    venous_pvs_diffusion = config["venous_pvs_diffusion"]
    parenchyma_diffusion = config["parenchyma_diffusion"]
    pvs_parenchyma_permability = config["pvs_parenchyma_permability"]
    pvs_csf_permability = config["pvs_csf_permability"]

    vein, vein_radii, vein_roots = read_vtk_network("mesh/networks/venes_smooth.vtk")
    vein_radii = as_P0_function(vein_radii)
    perm_vein  = 2*np.pi*vein_radii 
    area_vein  = np.pi *(vein_radii**2) 
    
    artery, artery_radii, artery_roots = read_vtk_network("mesh/networks/arteries_smooth.vtk")
    artery_radii = as_P0_function(artery_radii)
    perm_artery  = 2*np.pi*artery_radii 
    area_artery  = np.pi *(artery_radii**2)

    sas = Mesh()
    with XDMFFile('mesh/volmesh/mesh.xdmf') as f:
        f.read(sas)
        vol_subdomains = MeshFunction('size_t', sas, 3, 0)
        f.read(vol_subdomains, 'label')

    assert np.allclose(np.unique(vol_subdomains.array()), [1,2])

    # scale from mm to m
    [m.scale(1e-3) for m in [sas, vein, artery]]
    vein_radii.vector()[:] *= 1e-3
    artery_radii.vector()[:] *= 1e-3

    
    inlet = CompiledSubDomain("on_boundary && " + 
                            "+ (x[0] - m0)*(x[0] - m0) " + 
                            "+ (x[1] - m1)*(x[1] - m1) " + 
                            "+ (x[2] - m2)*(x[2] - m2) < r*r",
                            m0=inlet_midpoint[0], m1=inlet_midpoint[1],
                            m2=inlet_midpoint[2], r=inlet_radius)
    
    bm = MeshFunction("size_t", sas, 2, 0)
    inlet.mark(bm, 1)
    
    phi = pcws_constant(vol_subdomains, {1: Constant(1),  # csf
                                         2: Constant(ecs_share) # parenchyma
                                        })
    Ds = pcws_constant(vol_subdomains, {1: Constant(sas_diffusion),  # csf
                                        2: Constant(parenchyma_diffusion) # parenchyma
                                        })
    xi = pcws_constant(vol_subdomains, {1: Constant(pvs_csf_permability),  # csf
                                        2: Constant(pvs_parenchyma_permability) # parenchyma
                                        })
    xi.set_allow_extrapolation(True)
    Da = Constant(arterial_pvs_diffusion) 
    Dv = Constant(venous_pvs_diffusion)

    velocity_a = Constant(0.0)
    velocity_v = Constant(0.0)

    fa = Constant(0.0) 
    fv = Constant(0.0)

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
    artery_shape = xii.Circle(radius=artery_radii, degree=20,)
    ua, va = (xii.Average(x, artery, artery_shape) for x in (u, v))


    ##
    S3d = FunctionSpace(sas, "DG", 0)
    S1d = FunctionSpace(artery, "DG", 0)
    foo = interpolate(Constant(1), S3d)
    bar = TestFunction(S1d)
    
    vec = xii.ii_assemble((1/CellVolume(artery))*inner(xii.Average(foo, artery, artery_shape), bar)*dx_a)
    arterial_mask = Function(S1d)
    coefs = vec.get_local(); coefs[coefs < 0.99] = 0.0; vec.set_local(coefs)
    arterial_mask.vector()[:] = vec 
                                                                        

    dx_v = Measure('dx', domain=vein)
    vein_shape = xii.Circle(radius=vein_radii, degree=20)
    uv, vv = (xii.Average(x, vein, vein_shape) for x in (u, v)) 
    ##
    S3d = FunctionSpace(sas, "DG", 0)
    S1d = FunctionSpace(vein, "DG", 0)
    foo = interpolate(Constant(1), S3d)
    bar = TestFunction(S1d)
    
    vec = xii.ii_assemble((1/CellVolume(vein))*inner(xii.Average(foo, vein, vein_shape), bar)*dx_v)
    coefs = vec.get_local(); coefs[coefs < 0.99] = 0.0; vec.set_local(coefs)

    venous_mask = Function(S1d)
    venous_mask.vector()[:] = vec

    phia = xii.Average(phi, artery, artery_shape)
    phiv = xii.Average(phi, vein,vein_shape)
    xia = xii.Average(xi, artery, artery_shape)
    xiv = xii.Average(xi, vein, vein_shape)

    a = xii.block_form(W, 2)
    xia_int = interpolate(xi, Qa)
    xiv_int = interpolate(xi, Qv)

    a[0][0] = phi*(1/dt)*inner(u,v)*dx + phi*Ds*inner(grad(u), grad(v))*dx \
            + phia*xia*arterial_mask*(2*pi*perm_artery)*inner(ua, va)*dx_a \
            + phiv*xiv*venous_mask*(2*pi*perm_vein)*inner(uv, vv)*dx_v

    a[0][1] = -phia*xia*arterial_mask*(2*pi*perm_artery)*inner(pa, va)*dx_a
    a[0][2] = -phiv*xiv*venous_mask*(2*pi*perm_vein)*inner(pv, vv)*dx_v

    a[1][0] =-phia*xia*arterial_mask*(2*pi*perm_artery)*inner(qa, ua)*dx_a
    a[1][1] = (1/dt)*area_artery*inner(pa,qa)*dx + Da*area_artery*inner(grad(pa), grad(qa))*dx \
            - area_artery*inner(pa, dot(velocity_a,grad(qa)[0]))*dx  \
            + xia_int*arterial_mask*(2*pi*perm_artery)*inner(pa, qa)*dx

    a[2][0] = -phiv*xiv*venous_mask*(2*pi*perm_vein)*inner(qv, uv)*dx_v
    a[2][2] = (1/dt)*area_vein*inner(pv,qv)*dx + Dv*area_vein*inner(grad(pv), grad(qv))*dx \
            - area_vein*inner(pv, dot(velocity_v,grad(qv)[0]))*dx \
            + xiv_int*venous_mask*(2*pi*perm_vein)*inner(pv, qv)*dx 

    L = xii.block_form(W, 1)
    L[0]  = phi*(1/dt)*inner(u_i,v)*dx 
    
    L[1]  = (1/dt)*area_artery*inner(pa_i, qa)*dx + area_artery*inner(fa,qa)*dx
    L[2]  = (1/dt)*area_vein*inner(pv_i, qv)*dx + area_vein*inner(fv,qv)*dx


    AA, bb = map(xii.ii_assemble, (a, L))

    if config["bc_sas"] != "None":
        V_bcs  = [DirichletBC(V, config["bc_sas"], inlet)]
    else:
        V_bcs = []
    if config["bc_arteries"] != "None":
        Qa_bcs = [DirichletBC(Qa, config["bc_arteries"], inlet)]
    else:
        Qa_bcs = []
    if config["bc_venes"] != "None":
        Qv_bcs = [DirichletBC(Qa, config["bc_venes"], inlet)]
    else:
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
    #opts.setValue('ksp_rtol', 1E-40)
    opts.setValue('ksp_atol', 1E-12)   # |AX-b| < 1E-
    opts.setValue('pc_type', 'hypre')
    opts.setValue('ksp_initial_guess_nonzero', 1)

    ksp = PETSc.KSP().create()
    ksp.setOperators(A_, A_)
    ksp.setFromOptions()
    print('Start solve')
    t = 0.0 
    u_i.rename("c_sas", "time")
    pa_i.rename("c_artery", "time")
    pv_i.rename("c_vein", "time")

    results_dir = f"results/{modelname}/"
    os.makedirs(results_dir, exist_ok=True)
    pvdsas = File(results_dir + f'{modelname}_sas.pvd') 
    pvdarteries = File(results_dir + f'{modelname}_arteries.pvd') 
    pvdvenes = File(results_dir + f'{modelname}_venes.pvd') 
    files = (pvdsas, pvdarteries, pvdvenes)

    def write(sols, files, t):
        for s,f in zip(sols, files):
            f.write(s, t)

    write((u_i, pa_i, pv_i), files, 0.0)
    write((xia_int, xiv_int), (pvdarteries, pvdvenes), 0.0)


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
        wh[0].rename("c_sas", "time")
        wh[1].rename("c_artery", "time")
        wh[2].rename("c_vein", "time")
        write(wh, files, float(t))


if __name__ == "__main__":
    typer.run(run_simulation)
