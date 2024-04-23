from petsc4py import PETSc
from dolfin import *
import xii 
from solver import * 
import os
import numpy as np
import typer
from pathlib import Path
from plotting_utils import read_config


def write(sols, pvdfiles, t, hdffile=None):
    for s,f in zip(sols, pvdfiles):
        f.write(s, t)
    if hdffile is not None:
        for s in sols:
            hdffile.write(s, s.name(), t)

def run_simulation(configfile: str):

    config = read_config(configfile)
    modelname = Path(configfile).stem
    meshname = config["mesh"]

    results_dir = f"results/{modelname}/"
    os.makedirs(results_dir, exist_ok=True)

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
    pvs_basal_cistern_permability = config.get("pvs_basal_cistern_permability", pvs_csf_permability)
    pvs_ratio_venes = config["pvs_ratio_venes"]
    pvs_ratio_artery = config["pvs_ratio_artery"]
    beta = config["molecular_outflow_resistance"]

    vein, vein_radii, vein_roots = read_vtk_network("mesh/networks/venes_smooth.vtk")
    vein_radii = as_P0_function(vein_radii)
    vein_radii.vector()[:] *= pvs_ratio_venes
    perm_vein  = 2*np.pi*vein_radii
    area_vein  = np.pi*(vein_radii**2- (vein_radii/pvs_ratio_venes)**2)
    
    artery, artery_radii, artery_roots = read_vtk_network("mesh/networks/arteries_smooth.vtk")
    artery_radii = as_P0_function(artery_radii)
    artery_radii.vector()[:] *= pvs_ratio_artery
    perm_artery  = 2*np.pi*artery_radii 
    area_artery  = np.pi*(artery_radii**2 - (artery_radii/pvs_ratio_artery)**2)

    sas, vol_subdomains, sas_components = get_sas(meshname) 
    gdim = sas.geometric_dimension()
    
    assert np.allclose(np.unique(vol_subdomains.array()), [1,2])


    artery_shape = xii.Circle(radius=artery_radii, degree=20,)
    vein_shape = xii.Circle(radius=vein_radii, degree=20)
    artmarker = volmarker_to_networkmarker(vol_subdomains, artery, artery_shape)
    veinmarker = volmarker_to_networkmarker(vol_subdomains, vein, vein_shape)
    artmarker.rename("marker", "time")
    veinmarker.rename("marker", "time")
    vol_subdomains.rename("marker", "time")
    inlet_id = 2
    efflux_id = 3
    inlet = CompiledSubDomain("on_boundary && " + 
                            "+ (x[0] - m0)*(x[0] - m0) " + 
                            "+ (x[1] - m1)*(x[1] - m1) " + 
                            "+ (x[2] - m2)*(x[2] - m2) < r*r",
                            m0=inlet_midpoint[0], m1=inlet_midpoint[1],
                            m2=inlet_midpoint[2], r=inlet_radius)

    bm = MeshFunction("size_t", sas, 2, 0)
    inlet.mark(bm, inlet_id)
    efflux = CompiledSubDomain("on_boundary && x[2] > 0.2")
    efflux.mark(bm, efflux_id)

    basal_cistern = CompiledSubDomain("(x[0] - m0)*(x[0] - m0) " + 
                            "+ (x[1] - m1)*(x[1] - m1) " + 
                            "+ (x[2] - m2)*(x[2] - m2) < r*r",
                            m0=inlet_midpoint[0], m1=inlet_midpoint[1],
                            m2=inlet_midpoint[2], r=inlet_radius*1.5)
    cisternmarker = MeshFunction("size_t", artery, 1, 0)
    basal_cistern.mark(cisternmarker, 1)
    artmarker.array()[:] =  np.where(cisternmarker.array() * (artmarker.array()>0), 
                                     3, artmarker.array()[:])

    xi_dict = {0:Constant(0), 1: Constant(pvs_csf_permability), 
               2: Constant(pvs_parenchyma_permability),
               3: Constant(pvs_basal_cistern_permability)}
    phi_dict = {1: Constant(1),  2: Constant(ecs_share)}
    phi = pcws_constant(vol_subdomains, phi_dict)
    Ds = pcws_constant(vol_subdomains, {1: Constant(sas_diffusion),  # csf
                                        2: Constant(parenchyma_diffusion) # parenchyma
                                        })
    xi_a = pcws_constant(artmarker, xi_dict)
    xi_v = pcws_constant(veinmarker, xi_dict)

    Da = Constant(arterial_pvs_diffusion) 
    Dv = Constant(venous_pvs_diffusion)

    if "arterial_velocity_file" in config.keys():
        vel_file = config["arterial_velocity_file"]

        with XDMFFile(vel_file) as file:
            CG = FunctionSpace(artery, "CG", 1)
            velocity_a = Function(CG)
            file.read_checkpoint(velocity_a, "velocity")
        File(results_dir + f'{modelname}_flux.pvd') << velocity_a
    else:
        velocity_a = Constant(0)

    if "csf_velocity_file" in config.keys():
        vel_file = config["csf_velocity_file"]

        with XDMFFile(vel_file) as file:
            CG3 = VectorFunctionSpace(sas, "CG", 3)
            velocity_csf = Function(CG3)
            file.read_checkpoint(velocity_csf, "velocity")

    else:
        velocity_csf = Constant([0]*gdim)


    velocity_v = Constant([0]*gdim)

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
    ua, va = (xii.Average(x, artery, artery_shape) for x in (u, v))
                                               
    dx_v = Measure('dx', domain=vein)
    uv, vv = (xii.Average(x, vein, vein_shape) for x in (u, v)) 
    ds = Measure("ds", sas, subdomain_data=bm)
    # tangent vector
    tau = xii.TangentCurve(artery)
    a = xii.block_form(W, 2)

    a[0][0] = phi*(1/dt)*inner(u,v)*dx + phi*Ds*inner(grad(u), grad(v))*dx \
            + xi_a*(perm_artery)*inner(ua, va)*dx_a \
            + xi_v*(perm_vein)*inner(uv, vv)*dx_v \
            - inner(u, dot(velocity_csf, grad(v)))*dx \
            + beta*u*v*ds(efflux_id)

    a[0][1] = -xi_a*(perm_artery)*inner(pa, va)*dx_a
    a[0][2] = -xi_v*(perm_vein)*inner(pv, vv)*dx_v

    a[1][0] =-xi_a*(perm_artery)*inner(qa, ua)*dx_a
    a[1][1] = (1/dt)*area_artery*inner(pa,qa)*dx + Da*area_artery*inner(grad(pa), grad(qa))*dx \
            - area_artery*inner(pa, dot(velocity_a*tau,grad(qa)))*dx  \
            + xi_a*(perm_artery)*inner(pa, qa)*dx

    a[2][0] = -xi_v*(perm_vein)*inner(qv, uv)*dx_v
    a[2][2] = (1/dt)*area_vein*inner(pv,qv)*dx + Dv*area_vein*inner(grad(pv), grad(qv))*dx \
            - area_vein*inner(pv, dot(velocity_v,grad(qv)))*dx \
            + xi_v*(perm_vein)*inner(pv, qv)*dx 

    L = xii.block_form(W, 1)
    L[0]  = phi*(1/dt)*inner(u_i,v)*dx
    
    L[1]  = (1/dt)*area_artery*inner(pa_i, qa)*dx + area_artery*inner(fa,qa)*dx
    L[2]  = (1/dt)*area_vein*inner(pv_i, qv)*dx + area_vein*inner(fv,qv)*dx

    W_bcs = [[], [], []]
    expressions = []
    ds_a = Measure("ds", artery, subdomain_data=artery_roots)
    ds_v = Measure("ds", vein, subdomain_data=vein_roots)
    ds_i = [ds, ds_a, ds_v]
    dbcflag = False
    for i, dom in enumerate(["sas", "arteries", "venes"]):
        try:
            bc_conf = config["boundary_conditions"][dom]
        except KeyError:
            print(f"{dom} inlet boundary not configured, continuing...")
            continue
        expr_params = bc_conf.get("params", {})
        area = assemble(1*ds_i[i](inlet_id))
        expr = Expression(bc_conf["expr"], t=0, A=area, **expr_params, degree=1)
        expressions.append(expr)
        if bc_conf["type"] == "Dirichlet":
            W_bcs[i].append(DirichletBC(W[i], expr, ds_i[i].subdomain_data(), inlet_id))
            dbcflag = True
        if bc_conf["type"] == "Neumann":
            v = TestFunction(W[i])
            L[i] += expr*v*ds_i[i](inlet_id)

    AA, bb = map(xii.ii_assemble, (a, L))
    if dbcflag:
        AA, _, bc_apply_b = xii.apply_bc(AA, bb, bcs=W_bcs, return_apply_b=True)

    A_ = ksp_mat(xii.ii_convert(AA))

    opts = PETSc.Options() 
    opts.setValue('ksp_type', 'preonly')    
    #opts.setValue('ksp_view', None)
    #opts.setValue('ksp_view_eigenvalues', None)
    #opts.setValue('ksp_converged_reason', None)
    #opts.setValue('ksp_norm_type', 'unpreconditioned')
    #opts.setValue('ksp_monitor_true_residual', None)
    #opts.setValue('ksp_rtol', 1E-40)
    #opts.setValue('ksp_atol', 1E-12)   # |AX-b| < 1E-
    opts.setValue('pc_type', 'lu')
    opts.setValue("pc_factor_mat_solver_type", "mumps")
    opts.setValue("mat_mumps_icntl_4", "3")

    #opts.setValue('ksp_initial_guess_nonzero', 1)

    ksp = PETSc.KSP().create()
    ksp.setOperators(A_, A_)
    ksp.setFromOptions()
    print('Start solve')
    t = 0.0 
    u_i.rename("c_sas", "time")
    pa_i.rename("c_artery", "time")
    pv_i.rename("c_vein", "time")
    xi_a.rename("xi", "time")
    xi_v.rename("xi", "time")

    pvdsas = File(results_dir + f'{modelname}_sas.pvd') 
    pvdarteries = File(results_dir + f'{modelname}_arteries.pvd') 
    pvdvenes = File(results_dir + f'{modelname}_venes.pvd') 
    pvdfiles = (pvdsas, pvdarteries, pvdvenes)
    hdffile = HDF5File(sas.mpi_comm(), results_dir + f"{modelname}.hdf", "w")

    #write((u_i, pa_i, pv_i), pvdfiles, 0.0, hdffile=hdffile)
    #write((vol_subdomains, artmarker, veinmarker), pvdfiles, 0.0)
    #write((xi_a, xi_v), (pvdarteries, pvdvenes), 0.0)
    #write([phi], [pvdsas], 0.0)

    wh = xii.ii_Function(W)
    x_ = A_.createVecLeft()
    i = 0
    while t < T: 
        i += 1
        t += dt 
        for expr in expressions:
            expr.t = t
        print("time", t)
        bb = xii.ii_assemble(L)
        if dbcflag:
            bb = bc_apply_b(bb)
        b_ = ksp_vec(xii.ii_convert(bb))
        ksp.solve(b_, x_)
        # NOTE: solve(b_, ksp_vec(wh.vector())) segfault most likely because
        # of type incompatibility seq is expected and we have nest
        wh.vector()[:] = PETScVector(x_)
        u_i.assign(wh[0]) 
        pa_i.assign(wh[1]) 
        pv_i.assign(wh[2])

        wh[0].rename("c_sas", "time")
        wh[1].rename("c_artery", "time")
        wh[2].rename("c_vein", "time")
        #if i%config["output_frequency"] == 0:
           # write(wh, pvdfiles, float(t), hdffile)

if __name__ == "__main__":
    typer.run(run_simulation)
