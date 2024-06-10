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

def directsolve(a, L, bcs, a_prec, W):
    wh = Function(W)
    solve(a==L, wh, bcs=bcs, solver_parameters={"linear_solver":"mumps"})
    return wh

def iterativesolve(a, L, bcs, a_prec, W):
    A, b = assemble_system(a, L, bcs)
    B, _ = assemble_system(a_prec, L, bcs) 

    ns = Function(W)    # W is mixed FS
    nullspace_basis = ns.vector().copy()
    W.sub(1).dofmap().set(nullspace_basis, 1.0)
    nullspace_basis.apply('insert')
    nullspace = VectorSpaceBasis([nullspace_basis])
    nullspace.orthonormalize()

    as_backend_type(A).set_nullspace(nullspace)
    wh = Function(W)
    solver = PETScKrylovSolver()
    solver.set_operators(A, B)

    ksp = solver.ksp()
    opts = PETSc.Options()
    opts.setValue('ksp_type', 'minres')
    opts.setValue('ksp_minres_qlp', 1)
    opts.setValue('ksp_rtol', 1E-5)                
    opts.setValue('ksp_view_pre', None)
    opts.setValue('ksp_monitor_true_residual', None)                
    opts.setValue('ksp_converged_reason', None)
    # Specify that the preconditioner is block diagonal and customize
    # inverses of individual blocks
    opts.setValue('fieldsplit_0_ksp_type', 'preonly')
    opts.setValue('fieldsplit_0_pc_type', 'lu')
    opts.setValue('fieldsplit_0_pc_factor_mat_solver_type', 'mumps')
    opts.setValue('fieldsplit_1_ksp_type', 'preonly')
    opts.setValue('fieldsplit_1_pc_type', 'lu')
    opts.setValue('fieldsplit_1_pc_factor_mat_solver_type', 'mumps')
    pc = ksp.getPC()
    pc.setType(PETSc.PC.Type.FIELDSPLIT)
    is_V = PETSc.IS().createGeneral(W.sub(0).dofmap().dofs())
    is_Q = PETSc.IS().createGeneral(W.sub(1).dofmap().dofs())
    pc.setFieldSplitIS(('0', is_V), ('1', is_Q))
    pc.setFieldSplitType(PETSc.PC.CompositeType.ADDITIVE) 
    ksp.setUp()
    subksps = pc.getFieldSplitSubKSP()
    A0, B0 = subksps[0].getOperators()
    gdim = W.mesh().geometry().dim()
    A0.setBlockSize(gdim)
    B0.setBlockSize(gdim)
    pc.setFromOptions()
    ksp.setFromOptions()
    niters = solver.solve(wh.vector(), b)
    return wh

PETScOptions.set("mat_mumps_icntl_4", 3)  # mumps verbosity
#PETScOptions.set("mat_mumps_icntl_24", 1)  # null pivot detection
PETScOptions.set("mat_mumps_icntl_35", 1)  # BLR feature
PETScOptions.set("mat_mumps_cntl_7", 1e-8)  # BLR eps
#PETScOptions.set("mat_mumps_icntl_32", 1)  # forward elimination during solve (useful, but not passed on by petsc)
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

def map_on_global(uh, parentmesh, eps_digits=16):
    # map the solution back on the whole domain
    el = uh.function_space().ufl_element()
    V = uh.function_space()
    eldim = 1 if len(uh.ufl_shape)==0 else uh.ufl_shape[0]
    V_glob   = FunctionSpace(parentmesh, el)
    c_coords = np.round(V.tabulate_dof_coordinates()[::eldim,:], eps_digits)
    p_coords = np.round(V_glob.tabulate_dof_coordinates()[::eldim,:], eps_digits)
    idxmap = npi.indices(p_coords, c_coords, axis=0)
    uh_global = Function(V_glob)
    for i in range(eldim):
        uh_global.vector()[idxmap*eldim + i] = uh.vector()[i::eldim]
    return uh_global

def map_on_global2(uh, parentmesh, eps_digits=12):
    # map the solution back on the whole domain
    el = uh.function_space().ufl_element()
    V = uh.function_space()
    eldim = 1 if len(uh.ufl_shape)==0 else uh.ufl_shape[0]
    V_glob = FunctionSpace(parentmesh, el)
    
    # Increase precision for matching coordinates
    c_coords = np.round(V.tabulate_dof_coordinates()[::eldim, :], eps_digits)
    p_coords = np.round(V_glob.tabulate_dof_coordinates()[::eldim, :], eps_digits)
    
    try:
        idxmap = npi.indices(p_coords, c_coords, axis=0)
    except KeyError as e:
        # Handle missing coordinates
        idxmap = npi.indices(p_coords, c_coords, axis=0, missing='mask')
        if np.ma.is_masked(idxmap):
            missing_indices = np.where(np.ma.getmask(idxmap))[0]
            print(f"Some local coordinates are not present in the global coordinates. Missing indices: {missing_indices} length: {len(missing_indices)}, {len(c_coords)}")
            # Optionally, handle missing indices here
            idxmap = np.ma.filled(idxmap, fill_value= 0.0)  # or any other way to handle missing points

    uh_global = Function(V_glob)
    for i in range(eldim):
        if np.any(idxmap < 0):
            print(f"Warning: Some coordinates were not mapped correctly. Skipping these indices.")
        valid_idx = idxmap >= 0
        uh_global.vector()[idxmap[valid_idx] * eldim + i] = uh.vector()[i::eldim][valid_idx]

    return uh_global

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


def compute_sas_flow(configfile : str, elm:str = 'cr'):

    config = read_config(configfile)
    modelname = Path(configfile).stem
    meshname = config["mesh"]

    results_dir = f"results/csf_flow/{modelname}/"
    os.makedirs(results_dir, exist_ok=True)
    # get mesh 
    sas = Mesh()
    with XDMFFile(f'mesh/{meshname}/volmesh/mesh.xdmf') as f:
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
    
    # translate markers to the sas outer mesh 
    boundary_markers_outer = MeshFunction("size_t", sas_outer, sas_outer.topology().dim() - 1, 0)
    DomainBoundary().mark(boundary_markers_outer, NO_SLIP_ID)
    sas_outer.translate_markers(boundary_markers, (EFFLUX_ID, LV_INTERF_ID), marker_f=boundary_markers_outer)

    spinal_outlet = CompiledSubDomain("on_boundary && x[2] < zmin + eps",
                                       zmin=sas.coordinates()[:,2].min(), eps=1e-3)
    spinal_outlet.mark(boundary_markers_outer, SPINAL_OUTLET_IT)

    File("fpp.pvd") << boundary_markers_outer

    ds = Measure("ds", domain=sas_outer, subdomain_data=boundary_markers_outer)

    # Define function spaces for velocity and pressure
    cell = sas_outer.ufl_cell()  
    if elm == "BDM":
        Velm = FiniteElement('BDM', cell, 1)
        Qelm = FiniteElement('Discontinuous Lagrange', cell, 0) 
    elif elm == "cr":
        Velm = VectorElement('Crouzeix-Raviart', cell, 1)
        Qelm = FiniteElement('Discontinuous Lagrange', cell, 0)
    else: 
        Velm = VectorElement('Lagrange', cell, 3)
        Qelm = FiniteElement('Lagrange', cell, 2) 
    

    Q = FunctionSpace(sas_outer, Qelm)
    W = FunctionSpace(sas_outer, Velm * Qelm)

    u , p = TrialFunctions(W) 
    v , q = TestFunctions(W)
    n = FacetNormal(sas_outer)

    mu = Constant(config["mu"]) # units need to be checked 
    R = Constant(config["R"]) # 1e-5 Pa/(mm s)
    f = Constant([0]*gdim)

    LV_surface_area = assemble(1*ds(LV_INTERF_ID))
    total_production = config["production_rate"]
    g = total_production / LV_surface_area
    g_L_per_day = 1e3 *(60*60*24) * assemble(g*ds(LV_INTERF_ID))
    print(f"production rate: {g_L_per_day} L/day")

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

    if elm == 'cr':
        a += 2*(mu/avg(hF))*inner(jump(u), jump(v))*dS
    elif elm == 'BDM':
        a += Stabilization(sas_outer, u,v, mu, penalty)

    L = inner(f, v)*dx

    CG2 = VectorFunctionSpace(sas_outer, "CG", 2)
    inflow = get_normal_func(sas_outer)
    inflow.vector()[:] *= -g
    bcs = [DirichletBC(W.sub(0), Constant((0, 0, 0)), ds.subdomain_data(), NO_SLIP_ID),
           DirichletBC(W.sub(0), inflow, ds.subdomain_data(), LV_INTERF_ID)]
    if elm == 'BDM':
         a += -inner(dot(2*mu*D(v), n), Tangent(u, n))*ds(NO_SLIP_ID) \
            + 2*mu*(penalty/hF)*inner(Tangent(u, n), Tangent(v, n))*ds(NO_SLIP_ID) 
    if config["spinal_outflow_bc"] == "noslip":
        bcs += [DirichletBC(W.sub(0), Constant((0, 0, 0)), ds.subdomain_data(), SPINAL_OUTLET_IT)]
        if elm == 'BDM':
            a += -inner(dot(2*mu*D(v), n), Tangent(u, n))*ds(SPINAL_OUTLET_IT) \
                 + 2*mu*(penalty/hF)*inner(Tangent(u, n), Tangent(v, n))*ds(SPINAL_OUTLET_IT)  
    elif config["spinal_outflow_bc"] == "zeroneumann":
        pass
    else:
        raise Exception("spinal bc must be one of {noslip,zeroneumann}")

    a_prec = (inner(2*mu*grad(u), grad(v))*dx + inner(R*dot(u,n), dot(v,n))*ds(EFFLUX_ID)
                    + (1/mu)*inner(p, q)*dx) 

    #wh = iterativesolve(a, L, bcs, a_prec, W)
    wh = directsolve(a, L, bcs, a_prec, W)
    uh, ph = wh.split(deepcopy=True)[:]
    #assert np.isclose(assemble(inner(uh,-n)*ds(LV_INTERF_ID)), total_production, rtol=0.03)
    
    # project to CG2 and write for visualization
    uh2 = project(uh, CG2, solver_type="cg", preconditioner_type="hypre_amg")
    with XDMFFile(f'{results_dir}/csf_vis_v.xdmf') as xdmf:
        xdmf.write_checkpoint(uh2, "velocity")
    with XDMFFile(f'{results_dir}/csf_vis_p.xdmf') as xdmf:
        xdmf.write_checkpoint(ph, "pressure")
    
    from IPython import embed  
    embed()
    uh_global = map_on_global(uh, sas)
    dxglob = Measure("dx", sas, subdomain_data=label)

    assert np.isclose(assemble(inner(uh_global, uh_global)*dxglob(PARID)), 0)
    assert np.isclose(assemble(inner(uh_global, uh_global)*dxglob(CSFNOFLOWID)), 0)
    assert np.isclose(assemble(inner(uh, uh)*dx), assemble(inner(uh_global, uh_global)*dxglob))

    divu_global = assemble(div(uh_global)*div(uh_global)*dxglob)
    divu = assemble(div(uh)*div(uh)*dx)
    print(divu) 
    print(divu_global)
    #  assert np.isclose(divu, divu_global, rtol=0.2)

    with XDMFFile(f'{results_dir}/csf_v.xdmf') as xdmf:
        xdmf.write(sas)
        xdmf.write_checkpoint(uh_global, "velocity", append=True)

    with XDMFFile(f'{results_dir}/csf_p.xdmf') as xdmf:
        xdmf.write_checkpoint(ph, "pressure")

    umag = project(sqrt(inner(uh, uh)), FunctionSpace(sas_outer, "CG", 3),
                   solver_type="cg", preconditioner_type="hypre_amg")
    umean = assemble(umag*dx) / assemble(1*dx(domain=sas_outer))

    print(f"div u = {divu}")
    print(f"u max = {umag.vector().max()}")
    print(f"u umean = {umean}")

if __name__ == "__main__":
    typer.run(compute_sas_flow)