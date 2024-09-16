from dolfin import *
from solver import * 
import os
import numpy as np
from xii import *
from petsc4py import PETSc
import sympy as sp
import typer
import numpy_indexed as npi
from plotting_utils import read_config
from IPython import embed
from pathlib import Path
import yaml
from test_map_on_global_coords_shift import map_dg_on_global

def directsolve(a, L, bcs, a_prec, W):
    wh = Function(W)
    solve(a==L, wh, bcs=bcs, solver_parameters={"linear_solver":"mumps"})
    return wh

PETScOptions.set("mat_mumps_icntl_4", 3)  # mumps verbosity
#PETScOptions.set("mat_mumps_icntl_24", 1)  # null pivot detection
PETScOptions.set("mat_mumps_icntl_35", 1)  # BLR feature
PETScOptions.set("mat_mumps_cntl_7", 1e-8)  # BLR eps
#PETScOptions.set("mat_mumps_icntl_32", 1)  # forward elimination during solve (ptentially useful, but not passed on by petsc)
#PETScOptions.set("mat_mumps_icntl_22", 1)  # out-of-core to reduce memory
#PETScOptions.set("mat_mumps_icntl_11", 1)  # error analysis
#PETScOptions.set("mat_mumps_icntl_25", 2)  # turn on null space basis

# domain IDs
CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

# interface/boundary ids
EFFLUX_ID = 1
NO_SLIP_ID = 2
LV_INTERF_ID = 3
SPINAL_OUTLET_IT = 4
CSF_INTERF_ID = 5

def get_normal_func(mesh):
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(n, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L, "cg", "jacobi")
    return nh


def compute_sas_flow(configfile : str, elm:str = 'BDM'):

    config = read_config(configfile)
    modelname = Path(configfile).stem
    meshname = config["mesh"]

    order_k = 2

    results_dir = f"results/csf_flow/{modelname}/"
    os.makedirs(results_dir, exist_ok=True)
    # get mesh 
    sas = Mesh()
    with XDMFFile(meshname) as f:
        f.read(sas)
        gdim = sas.geometric_dimension()
        label = MeshFunction('size_t', sas, gdim, 0)
        f.read(label, 'label')

    sas_outer = EmbeddedMesh(label, [CSFID,LVID,V34ID]) 

    # create boundary markers 
    boundary_markers = MeshFunction("size_t", sas, sas.topology().dim() - 1)
    boundary_markers.set_all(0)

    # Sub domain for efflux route (mark whole boundary of the full domain) 
    efflux = CompiledSubDomain("on_boundary && x[2] > m", m=config["noslip_max_z"])
    efflux.mark(boundary_markers, EFFLUX_ID)
    
    # mark LV surface
    mark_internal_interface(sas, label, boundary_markers, LV_INTERF_ID,
                            doms=[LVID, PARID])
    
    mark_internal_interface(sas, label, boundary_markers, CSF_INTERF_ID,
                            doms=[CSFID, PARID])
    
    # translate markers to the sas outer mesh 
    boundary_markers_outer = MeshFunction("size_t", sas_outer, sas_outer.topology().dim() - 1, 0)
    DomainBoundary().mark(boundary_markers_outer, NO_SLIP_ID)
    sas_outer.translate_markers(boundary_markers,
                                (EFFLUX_ID, LV_INTERF_ID, CSF_INTERF_ID),
                                marker_f=boundary_markers_outer)

    spinal_outlet = CompiledSubDomain("on_boundary && x[2] < zmin + eps",
                                       zmin=sas.coordinates()[:,2].min(), eps=1e-3)
    spinal_outlet.mark(boundary_markers_outer, SPINAL_OUTLET_IT)

    File("fpp.pvd") << boundary_markers_outer

    ds = Measure("ds", domain=sas_outer, subdomain_data=boundary_markers_outer)

    # Define function spaces for velocity and pressure
    cell = sas_outer.ufl_cell()  
    Velm = FiniteElement('BDM', cell, order_k)
    Qelm = FiniteElement('Discontinuous Lagrange', cell, order_k - 1) 

    W = FunctionSpace(sas_outer, Velm * Qelm)

    u , p = TrialFunctions(W) 
    v , q = TestFunctions(W)
    n = FacetNormal(sas_outer)

    mu = Constant(config["mu"]) # units need to be checked 
    R = Constant(config["R"]) # 1e-5 Pa/(mm s)
    f = Constant([0]*gdim)

    LV_surface_area = assemble(1*ds(LV_INTERF_ID))
    PAR_surface_area = assemble(1*ds(CSF_INTERF_ID))
    SK_surface_area = assemble(1*ds(NO_SLIP_ID))
    g_LV = config["LV_inflow_rate"] / LV_surface_area
    g_par = config["tissue_inflow_rate"] / PAR_surface_area
    g_SK = config["skull_inflow_rate"] / SK_surface_area

    a = (inner(2*mu*sym(grad(u)), sym(grad(v)))*dx - inner(p, div(v))*dx
            -inner(q, div(u))*dx + inner(R*dot(u,n), dot(v,n))*ds(EFFLUX_ID)) 
   
    # NOTE: these are the jump operators from Krauss, Zikatonov paper.
    # Jump is just a difference and it preserves the rank 
    Jump = lambda arg: arg('+') - arg('-')
    # Average uses dot with normal and AGAIN MINUS; it reduces the rank
    Avg = lambda arg, n: Constant(0.5)*(dot(arg('+'), n('+')) - dot(arg('-'), n('-')))
    # Action of (1 - n x n)
    Tangent = lambda v, n: v - n*dot(v, n)    


    CellDiameter = CentroidDistance   # NOTE: also adjust the penalty parameter

    penalty = Constant(10.0)
    D = lambda v: sym(grad(v))

    def Stabilization(mesh, u, v, mu, penalty, consistent=True):
        '''Displacement/Flux Stabilization from Krauss et al paper'''
        n, hA = FacetNormal(mesh), avg(CellDiameter(mesh))
        

        if consistent:
            return (-inner(Avg(2*mu*D(u), n), Jump(Tangent(v, n)))*dS
                    -inner(Avg(2*mu*D(v), n), Jump(Tangent(u, n)))*dS
                    + 2*mu*(penalty/hA)*inner(Jump(Tangent(u, n)), Jump(Tangent(v, n)))*dS)
            
        # For preconditioning
        return 2*mu*(penalty/hA)*inner(Jump(Tangent(u, n)), Jump(Tangent(v, n)))*dS

    hF = CellDiameter(sas_outer)

    a += Stabilization(sas_outer, u,v, mu, penalty)

    L = inner(f, v)*dx
    LV_inflow = get_normal_func(sas_outer)
    LV_inflow.vector()[:] *= -g_LV
    PAR_inflow = get_normal_func(sas_outer)
    PAR_inflow.vector()[:] *= -g_par
    SK_inflow = get_normal_func(sas_outer)
    SK_inflow.vector()[:] *= -g_SK


    bcs = [
           DirichletBC(W.sub(0), SK_inflow, ds.subdomain_data(), NO_SLIP_ID),
           DirichletBC(W.sub(0), LV_inflow, ds.subdomain_data(), LV_INTERF_ID),
           DirichletBC(W.sub(0), PAR_inflow, ds.subdomain_data(), CSF_INTERF_ID)
           ]
    
    # add no-slip BC
    for bid in [NO_SLIP_ID, LV_INTERF_ID, CSF_INTERF_ID]:
        a += -inner(dot(2*mu*D(v), n), Tangent(u, n))*ds(bid) \
            - inner(dot(2*mu*D(u), n), Tangent(v, n))*ds(bid) \
            + 2*mu*(penalty/hF)*inner(Tangent(u, n), Tangent(v, n))*ds(bid)
        
    

    if config["spinal_outflow_bc"] == "noslip":
        bcs += [DirichletBC(W.sub(0), Constant((0, 0, 0)), ds.subdomain_data(), SPINAL_OUTLET_IT)]
        a += -inner(dot(2*mu*D(v), n), Tangent(u, n))*ds(SPINAL_OUTLET_IT) \
            - inner(dot(2*mu*D(u), n), Tangent(v, n))*ds(SPINAL_OUTLET_IT) \
                + 2*mu*(penalty/hF)*inner(Tangent(u, n), Tangent(v, n))*ds(SPINAL_OUTLET_IT)  
    elif config["spinal_outflow_bc"] == "zeroneumann":
        pass
    else:
        raise Exception("spinal bc must be one of {noslip,zeroneumann}")

    a_prec = (inner(2*mu*grad(u), grad(v))*dx + inner(R*dot(u,n), dot(v,n))*ds(EFFLUX_ID)
                    + (1/mu)*inner(p, q)*dx) 

    wh = directsolve(a, L, bcs, a_prec, W)
    uh, ph = wh.split(deepcopy=True)[:]
    
    uhdg = interpolate(uh, VectorFunctionSpace(sas_outer, "DG", order_k))
    assert np.isclose(assemble(div(uhdg)*div(uhdg)*dx), 0)
    with XDMFFile(f'{results_dir}/csf_vis_v.xdmf') as xdmf:
        xdmf.write_checkpoint(uhdg, "velocity")
    with XDMFFile(f'{results_dir}/csf_vis_p.xdmf') as xdmf:
        xdmf.write(sas_outer)
        xdmf.write_checkpoint(ph, "pressure", append=True)

    assert np.isclose(assemble(inner(uh,-n)*ds(LV_INTERF_ID)), config["LV_inflow_rate"], rtol=0.05)
    assert np.isclose(assemble(inner(uh,-n)*ds(CSF_INTERF_ID)), config["tissue_inflow_rate"], rtol=0.05)
    assert np.isclose(assemble(inner(uh,-n)*ds(NO_SLIP_ID)), config["skull_inflow_rate"], rtol=0.05)

    uh_global = Function(VectorFunctionSpace(sas, "DG", order_k))
    map_dg_on_global(uhdg, uh_global)
    dxglob = Measure("dx", sas, subdomain_data=label)

    #assert np.isclose(assemble(inner(uh_global, uh_global)*dxglob(PARID)), 0)

    divu_global_csf = assemble(div(uh_global)*div(uh_global)*dxglob(CSFID)) \
                    + assemble(div(uh_global)*div(uh_global)*dxglob(LVID)) \
                    + assemble(div(uh_global)*div(uh_global)*dxglob(V34ID))

    divu = assemble(sqrt(div(uh)*div(uh))*dx)

    assert np.isclose(divu, 0)

    with XDMFFile(f'{results_dir}/csf_v.xdmf') as xdmf:
        xdmf.write(sas)
        xdmf.write_checkpoint(uh_global, "velocity", append=True)

    ph_global_dg = Function(FunctionSpace(sas, "DG", order_k -1))
    map_dg_on_global(ph, ph_global_dg)
    with XDMFFile(f'{results_dir}/csf_p.xdmf') as xdmf:
        xdmf.write_checkpoint(ph_global_dg, "pressure")

    assert np.isclose(divu, divu_global_csf)
    #assert np.isclose(assemble(inner(uh_global, uh_global)*dxglob(CSFNOFLOWID)), 0)

    # collect key metrics:

    sas_vol = assemble(1*dx(domain=sas_outer))
    umean = assemble(sqrt(inner(uh, uh))*dx) / sas_vol
    umag = project(sqrt(inner(uh, uh)), FunctionSpace(sas_outer, "CG", order_k),
                   solver_type="cg", preconditioner_type="hypre_amg")
    umax = norm(umag.vector(), 'linf') 
    metrics = dict(umean=umean, umax=umax, divu=divu, 
                   pmax=ph.vector().max(), pmin=ph.vector().min())

    with open(f'{results_dir}/metrics.yml', 'w') as outfile:
        yaml.dump(metrics, outfile, default_flow_style=False)

if __name__ == "__main__":
    typer.run(compute_sas_flow)