from dolfin import *
from solver import * 
import os
import numpy as np
from xii import *
from petsc4py import PETSc
import typer
from plotting_utils import read_config
from IPython import embed
from pathlib import Path
import yaml

parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['optimize'] = True

# NOTE: these are the jump operators from Krauss, Zikatonov paper.
# Jump is just a difference and it preserves the rank 
Jump = lambda arg: arg('+') - arg('-')
# Average uses dot with normal and AGAIN MINUS; it reduces the rank
Avg = lambda arg, n: Constant(0.5)*(dot(arg('+'), n('+')) - dot(arg('-'), n('-')))
# Action of (1 - n x n)
Tangent = lambda v, n: v - n*dot(v, n)   
D = lambda v: sym(grad(v))

def Stabilization(mesh, u, v, mu, penalty, consistent=True):
    '''Displacement/Flux Stabilization from Krauss et al paper'''
    n, hA = FacetNormal(mesh), avg(CellDiameter(mesh))
    
    if consistent:
        return (-inner(Avg(mu*grad(u), n), Jump(Tangent(v, n)))*dS
                -inner(Avg(mu*grad(v), n), Jump(Tangent(u, n)))*dS
                + 2*mu*(penalty/hA)*inner(Jump(Tangent(u, n)), Jump(Tangent(v, n)))*dS)
        
    # For preconditioning
    return 2*mu*(penalty/hA)*inner(Jump(Tangent(u, n)), Jump(Tangent(v, n)))*dS


def directsolve(a, L, bcs, W):
    wh = Function(W)
    solve(a==L, wh, bcs=bcs, solver_parameters={"linear_solver":"mumps"})
    return wh

PETScOptions.set("mat_mumps_icntl_4", 3)  # mumps verbosity
#PETScOptions.set("mat_mumps_icntl_24", 1)  # null pivot detection
PETScOptions.set("mat_mumps_icntl_35", 1)  # BLR feature
#PETScOptions.set("mat_mumps_cntl_7", 1e-8)  # BLR eps
#PETScOptions.set("mat_mumps_icntl_22", 1)  # out-of-core to reduce memory
#PETScOptions.set("mat_mumps_icntl_25", 2)  # turn on null space basis#
PETScOptions.set("mat_mumps_icntl_28", 2)  # parallel ordering

# interface/boundary ids
UPPER_SKULL_ID = 1
LOWER_SKULL_ID = 2
LV_INTERF_ID = 3
PIA_ID = 4
SPINAL_CORD_ID = 5
SPINAL_OUTLET_ID = 6
CSF_NO_FLOW_CSF_ID = 7

def get_normal_func(mesh, scale=Constant(1)):
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "DG", 0)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(scale*n, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L, "cg", "jacobi")
    return nh

def get_inflow_bcs(W, ds, inflow_bcs):
    bcs = []
    for bids, infl in inflow_bcs:
        if infl == 0:
            infl_func = Constant((0,0,0))
        else:
            area = 0
            for bid in bids: area += assemble(1*ds(bid))
            infl_func = get_normal_func(W.mesh(), 
                scale=-Expression(infl, A=area, degree=3))
        for bid in bids:
            bcs += [DirichletBC(W.sub(0), infl_func, 
                                ds.subdomain_data(), bid)]
    return bcs

def get_BDM_problem(mesh, ds, order_k, mu, resistance_bcs, 
                    no_slip_bcs, inflow_bcs):
    
    cell = mesh.ufl_cell()  
    Velm = FiniteElement('BDM', cell, order_k)
    Qelm = FiniteElement('Discontinuous Lagrange',
                          cell, order_k  - 1) 

    W = FunctionSpace(mesh, Velm * Qelm)

    u , p = TrialFunctions(W) 
    v , q = TestFunctions(W)
    n = FacetNormal(mesh)
    penalty = Constant(10*order_k)
    CellDiameter = CentroidDistance
    hF = CellDiameter(mesh)
    mu = Constant(float(mu))
    f = Constant([0]*mesh.geometric_dimension())

    a = (inner(mu*grad(u), grad(v))*dx - inner(p, div(v))*dx
        - inner(q, div(u))*dx)
    
    a += Stabilization(mesh, u,v, mu, penalty)
    L = inner(f, v)*dx

    for bid, R in resistance_bcs:
        a += inner(Constant(float(R))*dot(u,n), dot(v,n))*ds(bid)

    for bid in no_slip_bcs:
        a += -inner(dot(mu*grad(v), n), Tangent(u, n))*ds(bid) \
            - inner(dot(mu*grad(u), n), Tangent(v, n))*ds(bid) \
            + 2*mu*(penalty/hF)*inner(Tangent(u, n), Tangent(v, n))*ds(bid)

    bcs = get_inflow_bcs(W, ds, inflow_bcs)
    
    return a, L, bcs, W

def get_TH_problem(mesh, ds, order_k, mu, resistance_bcs, 
                    no_slip_bcs, inflow_bcs):
    cell = mesh.ufl_cell()  
    Velm = VectorElement('CG', cell, order_k)
    Qelm = FiniteElement('CG', cell, order_k  - 1) 

    W = FunctionSpace(mesh, Velm * Qelm)

    u , p = TrialFunctions(W) 
    v , q = TestFunctions(W)
    n = FacetNormal(mesh)
    mu = Constant(mu)
    zero_vec = Constant([0]*mesh.geometric_dimension())

    a = (inner(mu*grad(u), grad(v))*dx - inner(p, div(v))*dx
        - inner(q, div(u))*dx)
    L = inner(zero_vec, v)*dx

    for bid, R in resistance_bcs:
        a += inner(R*dot(u,n), dot(v,n))*ds(bid)

    bcs = []
    flatten = lambda l: [x for el in l for x in el]
    inflow_bc_ids = flatten([bids for bids, infl in inflow_bcs])
    for bid in set(no_slip_bcs) - set(inflow_bc_ids):
        bcs += [DirichletBC(W.sub(0), zero_vec, 
                            ds.subdomain_data(), bid)]
    bcs += get_inflow_bcs(W, ds, inflow_bcs)
    return a, L, bcs, W

def compute_sas_flow(configfile : str):
    config = read_config(configfile)
    modelname = Path(configfile).stem
    meshname = config["mesh"]
    if config["discretization"] == "BDM": 
        parameters["ghost_mode"] = "shared_facet"
    order_k = config.get("order", 2)
    results_dir = f"results/csf_flow/{modelname}/"
    os.makedirs(results_dir, exist_ok=True)
    # get mesh 
    mesh = Mesh(MPI.comm_world)
    with XDMFFile(meshname) as f:
        f.read(mesh)
        sm = MeshFunction("size_t", mesh, 3, 0)
        f.read(sm, "f")
        gdim = mesh.geometric_dimension()

    bm = MeshFunction("size_t", mesh, 2, 0)
    with XDMFFile(meshname.replace(".xdmf","_facets.xdmf")) as f:
        f.read(bm, "f")
    bm.array()[bm.array()[:]==CSF_NO_FLOW_CSF_ID] = PIA_ID
    ds = Measure("ds", domain=mesh, subdomain_data=bm)

    assert assemble(1*ds(0)) == 0
    assert assemble(1*ds(CSF_NO_FLOW_CSF_ID)) == 0

    if config["discretization"] == "BDM":
        a, L, bcs, W = get_BDM_problem(mesh, ds, order_k, config["mu"],
                                       config["resistance_bcs"],
                                       config["no_slip_bcs"],
                                       config["inflow_bcs"])
    elif config["discretization"] == "TH":
        a, L, bcs, W = get_TH_problem(mesh, ds, order_k, config["mu"], 
                                      config["resistance_bcs"],
                                      config["no_slip_bcs"], 
                                      config["inflow_bcs"])
    else: print("choose BDM or TH discretization!"); exit()

    wh = directsolve(a, L, bcs, W)
    uh, ph = wh.split(deepcopy=True)[:]
    
    uhdg = interpolate(uh, VectorFunctionSpace(mesh, "DG", order_k))
    if ph.ufl_element().family() == 'Discontinuous Lagrange':
        assert np.isclose(assemble(div(uhdg)*div(uhdg)*dx), 0)
    else:
        ph = interpolate(ph, FunctionSpace(mesh, "DG", order_k - 1))

    with XDMFFile(f'{results_dir}/csf_v.xdmf') as xdmf:
        xdmf.write_checkpoint(uhdg, "velocity")
    with XDMFFile(f'{results_dir}/csf_p.xdmf') as xdmf:
        xdmf.write_checkpoint(ph, "pressure")

    with HDF5File(MPI.comm_world, f'{results_dir}/flow.hdf','w') as f:
        f.write(mesh, "mesh")
        f.write(ph, "pressure")
        f.write(uhdg, "velocity")
        f.write(as_P0_function(sm), "label")

    # collect key metrics:
    sas_vol = assemble(1*dx(domain=mesh))
    umean = assemble(sqrt(inner(uh, uh))*dx) / sas_vol
    umag = project(sqrt(inner(uh, uh)), FunctionSpace(mesh, "CG", order_k),
                   solver_type="cg", preconditioner_type="hypre_amg")
    umax = norm(umag.vector(), 'linf') 
    divu = assemble(sqrt(div(uhdg)*div(uhdg))*dx)
    metrics = dict(umean=umean, umax=umax, divu=divu, 
                   pmax=ph.vector().max(), pmin=ph.vector().min())

    if MPI.comm_world.rank == 0:
        with open(f'{results_dir}/metrics.yml', 'w') as outfile:
            yaml.dump(metrics, outfile, default_flow_style=False)


def pipetest():
    mesh = UnitSquareMesh(10, 10)
    walls = CompiledSubDomain("x[1] < DOLFIN_EPS || x[1] >= 1 - DOLFIN_EPS")
    inlet = CompiledSubDomain("x[0] < DOLFIN_EPS")
    bm = MeshFunction("size_t", mesh, 1, 0)
    inlet.mark(bm, 1)
    walls.mark(bm, 2)
    ds = Measure("ds", mesh, subdomain_data=bm)
    mu = 0.23
    order_k = 2
    a, L, bcs, W = get_BDM_problem(mesh, ds, order_k, mu, 
                                    [],
                                    [2], 
                                    [[[1], "-pow(2*x[1]- 1, 2) + 1"], [[2], "0"]])

    u,p = directsolve(a, L, bcs, W).split(deepcopy=True)
    u = interpolate(u, VectorFunctionSpace(mesh, "DG", order_k))
    from vtk_adapter import plot
    plot(u, "p.png")

if __name__ == "__main__":
    typer.run(compute_sas_flow)
    #mms()
