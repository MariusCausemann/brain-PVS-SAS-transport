from petsc4py import PETSc
from dolfin import *
import xii 
from solver import * 
import os
import numpy as np
import typer
from pathlib import Path
from plotting_utils import read_config
from test_map_on_global_coords_shift import map_dg_on_global, map_kdtree
from IPython import embed
import yaml

def write(sols, pvdfiles, t, hdffile=None):
    for s,f in zip(sols, pvdfiles):
        f.write(s, t)
    if hdffile is not None:
        for s in sols:
            hdffile.write(s, s.name(), t)

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

UPPER_SKULL_ID = 1
LOWER_SKULL_ID = 2
LV_INTERF_ID = 3
PIA_ID = 4
SPINAL_CORD_ID = 5
SPINAL_OUTLET_ID = 6
CSF_NO_FLOW_CSF_ID = 7

def get_subdomain_dofs(V, subdomains, subd_id):
    dm = V.dofmap()
    subd_dofs = np.unique(np.hstack(
        [dm.cell_dofs(c.index()) for c in SubsetIterator(subdomains, subd_id)]))
    return subd_dofs

def read_pvs_velocity(vel_file, netw):
    from evaluate_pvs_flow import sig_to_space
    print(f"reading arterial PVS velocity from {vel_file}")
    with HDF5File(MPI.comm_world, vel_file,'r') as f:
        sig = f.attributes("/u").to_dict()["signature"]
        vel = Function(sig_to_space(sig, netw))
        f.read(vel, "u")
    # make sure the velocity field is tangential to the network
    tau = xii.TangentCurve(netw)
    un = vel - tau*inner(vel, tau)
    assert assemble(inner(un,un)*dx) < 1e-16
    return vel 


def run_simulation(configfile: str):

    config = read_config(configfile)
    modelname = Path(configfile).stem
    meshname = config["mesh"]
    metrics = dict()
    results_dir = f"results/{modelname}/"
    os.makedirs(results_dir, exist_ok=True)

    dt = config["dt"]
    T = config["T"]

    ecs_share = config["ecs_share"]
    sas_diffusion = config["sas_diffusion"]
    arterial_pvs_diffusion = config["arterial_pvs_diffusion"]
    venous_pvs_diffusion = config["venous_pvs_diffusion"]
    parenchyma_diffusion = config["parenchyma_diffusion"]
    arterial_pvs_parenchyma_permability = config["arterial_pvs_parenchyma_permability"]
    arterial_pvs_csf_permability = config["arterial_pvs_csf_permability"]
    venous_pvs_parenchyma_permability = config["venous_pvs_parenchyma_permability"]
    venous_pvs_csf_permability = config["venous_pvs_csf_permability"]
    pvs_ratio_venes = config["pvs_ratio_venes"]
    pvs_ratio_artery = config["pvs_ratio_artery"]
    beta = config["molecular_outflow_resistance"]
    beta_csf_par = config["par_csf_permability"]
    root_xi_factor = config["root_pvs_permability_factor"]
    c_init = config["initial_concentration"]

    vein, vein_radii, vein_roots = read_vtk_network("mesh/networks/venes_smooth.vtk", rescale_mm2m=False)
    vein_radii = as_P0_function(vein_radii)
    vein_radii.vector()[:] *= pvs_ratio_venes
    perm_vein  = 2*np.pi*vein_radii
    area_vein  = np.pi*(vein_radii**2 - (vein_radii/pvs_ratio_venes)**2)
    
    artery, artery_radii, artery_roots = read_vtk_network("mesh/networks/arteries_smooth.vtk", rescale_mm2m=False)
    artery_radii = as_P0_function(artery_radii)
    artery_radii.vector()[:] *= pvs_ratio_artery
    perm_artery  = 2*np.pi*artery_radii 
    area_artery  = np.pi*(artery_radii**2 - (artery_radii/pvs_ratio_artery)**2)

    mesh, vol_subdomains = get_mesh(meshname) 
    bm = MeshFunction("size_t", mesh, 2, 0)
    with XDMFFile(meshname.replace(".xdmf","_facets.xdmf")) as f:
        f.read(bm, "f")
    bm.array()[bm.array()[:]==CSF_NO_FLOW_CSF_ID] = 0
    gdim = mesh.geometric_dimension()
    
    label = MeshFunction('size_t', mesh, gdim, 0)
    label.array()[:] = vol_subdomains.array()[:]
    assert np.allclose(np.unique(vol_subdomains.array()), [CSFID, PARID, LVID, V34ID, CSFNOFLOWID])
    vol_subdomains.array()[np.isin(vol_subdomains.array(), [LVID, V34ID, CSFNOFLOWID])] = CSFID
    File(results_dir + "subdomains.pvd") << vol_subdomains

    artery_shape = xii.Circle(radius=artery_radii, degree=40, quad_rule='midpoint')
    vein_shape = xii.Circle(radius=vein_radii, degree=40, quad_rule='midpoint')
    csf_par_weights = {PARID:0.2, CSFID:0.8}
    artmarker = volmarker_to_networkmarker(vol_subdomains, artery, artery_shape,
                                         filename=results_dir + "arttagshares.pvd",
                                         weights=csf_par_weights)
    veinmarker = volmarker_to_networkmarker(vol_subdomains, vein, vein_shape,
                                             filename=results_dir +  "ventagshares.pvd",
                                             weights=csf_par_weights)

    priority = xii.InterfaceResolution(subdomains=vol_subdomains,
                                       resolve_conflicts={(CSFID, PARID):
                                                           csf_par_weights[CSFID]/csf_par_weights[PARID]})
    artmarker.rename("marker", "time")
    veinmarker.rename("marker", "time")
    vol_subdomains.rename("marker", "time")

    art_xi_dict = {CSFID: arterial_pvs_csf_permability, 
               PARID: arterial_pvs_parenchyma_permability}
    ven_xi_dict = {CSFID: venous_pvs_csf_permability, 
               PARID: venous_pvs_parenchyma_permability}
    phi_dict = {CSFID: 1,  PARID: ecs_share}

    phi = pcws_constant(vol_subdomains, phi_dict)
    Ds = pcws_constant(vol_subdomains, {CSFID: sas_diffusion,  # csf
                                        PARID: parenchyma_diffusion # parenchyma
                                        })
    
    def get_xi(marker, xi_dict, roots):
        netw = marker.mesh()
        netw_xi_root =  Function(FunctionSpace(netw, "CG", 1))
        idx = map_kdtree(netw.coordinates(), FunctionSpace(netw, "CG", 1).tabulate_dof_coordinates())
        netw_xi_root.vector()[:] = np.where(roots.array()[:] > 0, root_xi_factor, 1)[idx]
        xi = pcws_constant(marker, xi_dict)
        xi.vector()[:] *= interpolate(netw_xi_root, xi.function_space()).vector()[:]
        return xi

    xi_a = get_xi(artmarker, art_xi_dict, artery_roots)
    xi_v = get_xi(veinmarker, ven_xi_dict, vein_roots)

    Da = Constant(arterial_pvs_diffusion) 
    Dv = Constant(venous_pvs_diffusion)

    velocity_a, velocity_v = Constant([0]*gdim), Constant([0]*gdim)

    for afile in config["arterial_velocity_file"]:
        velocity_a += read_pvs_velocity(afile, artery)
    for vfile in config["venous_velocity_file"]:
        velocity_v += read_pvs_velocity(vfile, vein)
        
    if "csf_velocity_file" in config.keys():
        print("reading CSF velocity from file...")
        vel_file = config["csf_velocity_file"]
        vel_mesh = Mesh()
        with HDF5File(MPI.comm_world, vel_file,'r') as f:
            f.read(vel_mesh, "mesh", False)
            v_elem = eval(f.attributes("/velocity").to_dict()["signature"])
            V = FunctionSpace(vel_mesh, v_elem)
            velocity_csf = Function(V)
            f.read(velocity_csf, "velocity")
            velocity_csf = map_dg_on_global(velocity_csf, parentmesh=mesh)
            velocity_csf.rename("v", "v")
            dx_s = Measure("dx", mesh, subdomain_data=vol_subdomains)
            assert assemble(sqrt(div(velocity_csf)*div(velocity_csf))*dx) < 1e-10
            assert assemble(inner(velocity_csf, velocity_csf)*dx_s(PARID)) < 1e-14
            with XDMFFile(results_dir + "csf_v.xdmf") as outfile:
                outfile.write_checkpoint(velocity_csf,"v", 0, append=False)
    else:
        velocity_csf = Constant([0]*gdim)


    if "csf_dispersion_file" in config.keys():
        dispfiles = config["csf_dispersion_file"]
        if isinstance(dispfiles, str): dispfiles = [dispfiles]
        dispweights = config.get("csf_dispersion_weights", [1]*len(dispfiles))

        sm = xii.EmbeddedMesh(label, [CSFID, LVID, V34ID])
        smDG0 = FunctionSpace(sm, "DG", 1)
        R = Function(smDG0)
        for dispf,w in zip(dispfiles, dispweights):
            print(f"reading CSF dispersion from file: {dispf}")
            with XDMFFile(dispf) as file:
                Ri = Function(smDG0)
                file.read_checkpoint(Ri, "R")
                R.vector()[:] += Ri.vector()[:] * w
        metrics["R_mean"] = assemble(R*dx) / assemble(1*dx(domain=sm))
        metrics["R_max"] = R.vector().max()
        metrics["R_min"] = R.vector().min()
        R = map_dg_on_global(R, parentmesh=mesh)
        Ds *= (1 + R)


    fa = Constant(0.0) 
    fv = Constant(0.0)

    V = FunctionSpace(mesh, 'DG', 1)
    Qa = FunctionSpace(artery, 'CG', 1)
    Qv = FunctionSpace(vein, 'CG', 1)
    W = [V, Qa, Qv]

    u, pa, pv = map(TrialFunction, W)
    v, qa, qv = map(TestFunction, W)
     
    # initial conditions 
    u_o  = Constant(c_init)
    pa_o = Constant(0.0) 
    pv_o = Constant(0.0) 

    u_i  = interpolate(u_o, V) 
    pa_i = interpolate(pa_o, Qa) 
    pv_i = interpolate(pv_o, Qv)
    # Things for restriction
    dx_a = Measure('dx', domain=artery)
    ua, va = (xii.Average(x, artery, artery_shape, normalize=True,
     resolve_interfaces=priority) for x in (u, v))
                                               
    dx_v = Measure('dx', domain=vein)
    uv, vv = (xii.Average(x, vein, vein_shape, normalize=True,
     resolve_interfaces=priority) for x in (u, v)) 

    ds = Measure("ds", mesh, subdomain_data=bm)
    # tangent vector
    a = xii.block_form(W, 2)

    eta = 1.0 
    alpha = Constant(1e3)

    dx_s = Measure("dx", mesh, subdomain_data=vol_subdomains)
    dS = Measure("dS", mesh, subdomain_data=bm)

    def a_dg_adv_diff(u,v) :

        # Bilinear form for DG advection diffusion
        n = FacetNormal(mesh)
        h = CellDiameter(mesh)
        b = velocity_csf
        dSi = dS(0)
        a_int = dot(grad(v), Ds*phi*grad(u))*dx - dot(grad(v),b*u)*dx_s(CSFID)

        wavg = lambda k: 2*k("+")*k("-") / (k("+") + k("-"))
        DF = wavg(Ds*phi)

        a_fac = (alpha/avg(h))*DF*dot(jump(u, n), jump(v, n))*dSi \
                - dot(avg(grad(u))*DF, jump(v, n))*dSi \
                - dot(jump(u, n), avg(grad(v)) * DF)*dSi \
                + beta_csf_par*jump(u)*jump(v)*dS(PIA_ID) \
                + beta_csf_par*jump(u)*jump(v)*dS(LV_INTERF_ID)
        
        a_vel = dot(dot(b("+"),n('+'))*avg(u), jump(v))*dSi \
            + (eta/2)*dot(abs(dot(b("+"),n('+')))*jump(u), jump(v))*dSi

        a = a_int + a_fac + a_vel
        return a

    def supg_stabilization(u, v, vel):
        
        mesh = u.function_space().mesh()
        hK = CellDiameter(mesh)

        beta = conditional(lt(sqrt(inner(vel, vel)), Constant(1E-10)),
                            Constant(0),
                            1/sqrt(inner(vel, vel)))

        dx5 = Measure("dx", metadata={'quadrature_degree': 8})

        return beta*hK*inner(dot(vel, grad(u)), dot(vel, grad(v)))*dx5

    a[0][0] = phi*(1/dt)*inner(u,v)*dx \
            + a_dg_adv_diff(u,v) \
            + beta*u*v*ds(UPPER_SKULL_ID) \
            + xi_a*(perm_artery)*inner(ua, va)*dx_a \
            + xi_v*(perm_vein)*inner(uv, vv)*dx_v

    a[0][1] = -xi_a*(perm_artery)*inner(pa, va)*dx_a
    a[0][2] = -xi_v*(perm_vein)*inner(pv, vv)*dx_v

    a[1][0] =-xi_a*(perm_artery)*inner(qa, ua)*dx_a
    a[1][1] = (1/dt)*area_artery*inner(pa,qa)*dx + Da*area_artery*inner(grad(pa), grad(qa))*dx \
            - area_artery*inner(pa, dot(velocity_a,grad(qa)))*dx  \
            + xi_a*(perm_artery)*inner(pa, qa)*dx \
            + supg_stabilization(pa, area_artery*qa, velocity_a) 


    a[2][0] = -xi_v*(perm_vein)*inner(qv, uv)*dx_v
    a[2][2] = (1/dt)*area_vein*inner(pv,qv)*dx + Dv*area_vein*inner(grad(pv), grad(qv))*dx \
            - area_vein*inner(pv, dot(velocity_v,grad(qv)))*dx \
            + xi_v*(perm_vein)*inner(pv, qv)*dx \
            + supg_stabilization(pv, area_vein*qv, velocity_v) 

    L     = xii.block_form(W, 1)
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
        area = assemble(1*ds_i[i](SPINAL_OUTLET_ID))
        expr = Expression(bc_conf["expr"], t=0, A=area, **expr_params, degree=1)
        expressions.append(expr)
        if bc_conf["type"] == "Dirichlet":
            W_bcs[i].append(DirichletBC(W[i], expr, ds_i[i].subdomain_data(), SPINAL_OUTLET_ID))
            dbcflag = True
        if bc_conf["type"] == "Neumann":
            v = TestFunction(W[i])
            L[i] += expr*v*ds_i[i](SPINAL_OUTLET_ID)

    AA, bb = map(xii.ii_assemble, (a, L))
    if dbcflag:
        AA, _, bc_apply_b = xii.apply_bc(AA, bb, bcs=W_bcs, return_apply_b=True)

   
    A_ = ksp_mat(xii.ii_convert(AA))

    opts = PETSc.Options() 
    opts.setValue('ksp_type', 'preonly')    
    opts.setValue('pc_type', 'lu')
    opts.setValue("pc_factor_mat_solver_type", "mumps")
    opts.setValue("mat_mumps_icntl_4", "3")
    opts.setValue("mat_mumps_icntl_35", 1)
    opts.setValue("mat_mumps_cntl_7",  1e-8)  # BLR eps

    ksp = PETSc.KSP().create()
    ksp.setOperators(A_, A_)
    ksp.setFromOptions()
    print('Start solve')
    t = 0.0 
    u_i.rename("c", "c")
    pa_i.rename("c", "c")
    pv_i.rename("c", "c")
    xi_a.rename("xi_a", "xi_a")

    xdmfsas = XDMFFile(results_dir + f'{modelname}_sas.xdmf') 
    xdmfart = XDMFFile(results_dir + f'{modelname}_artery.xdmf') 
    xdmfven = XDMFFile(results_dir + f'{modelname}_vein.xdmf') 
    xdmfsas.parameters["flush_output"] = True
    xdmfart.parameters["flush_output"] = True
    xdmfsas.parameters["flush_output"] = True

    pvdarteries = File(results_dir + f'{modelname}_artery.pvd') 
    pvdvenes = File(results_dir + f'{modelname}_vein.pvd') 
    pvdfiles = (pvdarteries, pvdvenes)

    pvdarteries.write(pa_i, t)
    pvdvenes.write(pv_i, t)
    xdmfsas.write_checkpoint(u_i, "c", t, append=False)
    xdmfart.write_checkpoint(pa_i, "c", t, append=False)
    xdmfven.write_checkpoint(pv_i, "c", t, append=False)

    #hdffile = HDF5File(mesh.mpi_comm(), results_dir + f"{modelname}.hdf", "w")

    #write((u_i, pa_i, pv_i), pvdfiles, 0.0, hdffile=hdffile)
    write((artmarker, veinmarker), pvdfiles, 0.0)
    write((xi_a, xi_v), pvdfiles, 0.0)
    #write([phi], [pvdsas], 0.0)

    wh = xii.ii_Function(W)
    x_ = A_.createVecLeft()
    i = 0
    par_dofs = get_subdomain_dofs(V, vol_subdomains, PARID)
    csf_dofs = get_subdomain_dofs(V, vol_subdomains, CSFID)

    for dom in ["par", "csf", "ven", "art"]:
        metrics[f"{dom}_min"] = [] 
        metrics[f"{dom}_max"] = [] 
        metrics[f"{dom}_mean"] = [] 

    while t < T: 
        i += 1
        t += dt 
        for expr in expressions:
            expr.t = t
        
        bb = xii.ii_assemble(L)
        if dbcflag:
            bb = bc_apply_b(bb)
        b_ = ksp_vec(xii.ii_convert(bb))
        ksp.solve(b_, x_)
        
        x2 = x_.norm(2)

        print("time", t, '|b|', b_.norm(2), '|x|', x2)

        # NOTE: solve(b_, ksp_vec(wh.vector())) segfault most likely because
        # of type incompatibility seq is expected and we have nest
        wh.vector()[:] = PETScVector(x_)
        u_i.assign(wh[0]) 
        pa_i.assign(wh[1]) 
        pv_i.assign(wh[2])

        wh[0].rename("c", "c")
        wh[1].rename("c", "c")
        wh[2].rename("c", "c")

        for dom, v, dofs, dxd in zip(["par", "csf", "art", "ven"],
                                [u_i, u_i, pa_i, pv_i],
                                [par_dofs, csf_dofs, None, None],
                                [dx_s(PARID), dx_s(CSFID), dx_a, dx_v]):
            metrics[f"{dom}_min"].append(v.vector().get_local()[dofs].min())
            metrics[f"{dom}_max"].append(v.vector().get_local()[dofs].max())
            mean = assemble(v*dxd) / assemble(1*dxd)
            metrics[f"{dom}_mean"].append(mean)

        if i%config["output_frequency"] == 0:
            #write(wh, pvdfiles, float(t), hdffile)
            xdmfsas.write_checkpoint(wh[0], "c", t, append=True)
            xdmfart.write_checkpoint(wh[1], "c", t, append=True)
            xdmfven.write_checkpoint(wh[2], "c", t, append=True)
            pvdarteries.write(pa_i, t)
            pvdvenes.write(pv_i, t)
        if np.isnan(x2) or x2 > 1e9:
            raise OverflowError(f"mumps produced nans at time {t}: |x| = {x2}")

    with open(results_dir + f'{modelname}_metrics.yml', 'w') as outfile:
        yaml.dump(metrics, outfile, default_flow_style=False)

if __name__ == "__main__":
    typer.run(run_simulation)
