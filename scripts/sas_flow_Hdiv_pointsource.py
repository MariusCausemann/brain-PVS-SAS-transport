from dolfin import *
from solver import * 
from xii import *
from xii.meshing.make_mesh_cpp import make_mesh
import os
import numpy as np
import typer
from plotting_utils import read_config
from pathlib import Path
import yaml
from pykdtree.kdtree import KDTree
from collections import namedtuple
# Represent flux out of the tube through a disk of radius passing through plane
# defined by x and normal.
FluxPointSource = namedtuple('FluxPointSource', ('x', 'normal', 'radius', 'value'))
m3s_to_mlday = 1e6*(60*60*24)

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


def point_source_mesh(sources, h):
    '''Auxiliary 1d mesh for Lagrange multipliers'''
    # NOTE: the mesh is 1d because of xii.Avarage
    coordinates, cells = [], []
    for (i, source) in enumerate(sources):
        v0 = source.x
        # Now we make a short cell
        n = source.normal
        n = n/np.linalg.norm(n)
        v1 = v0*h
        coordinates.extend((v0, v1))
        cells.append((2*i, 2*i+1))
    coordinates, cells = map(np.array, (coordinates, cells))
    print(coordinates)
    print(cells)
    mesh = make_mesh(coordinates, cells, 1, 3)
    x = mesh.coordinates()

    vertex_f = MeshFunction('size_t', mesh, 0, 0)
    for (k, source) in enumerate(sources, 1):
        dist = np.linalg.norm(x - source.x, 2, axis=1)
        i = np.argmin(dist)
        assert dist[i] < 1E-13
        vertex_f[i] = k

    return vertex_f


def load_pvs_sources(pvs_flow_file):
    from solver import read_vtk_network
    from evaluate_pvs_flow import sig_to_space

    mesh, artery_radii, artery_roots = read_vtk_network("mesh/networks/arteries_smooth.vtk", rescale_mm2m=False)

    with HDF5File(MPI.comm_world, pvs_flow_file,'r') as f:
        u_sig = f.attributes("/u").to_dict()["signature"]
        u = Function(sig_to_space(u_sig, mesh))
        f.read(u, "u")
    return u, artery_radii, artery_roots

def get_pvs_point_sources(sas_mesh_func):
    pvs_flow_file = "results/pvs_flow_peristaltic/vasomotion-strong/pvs_flow.hdf"
    u, artery_radii, artery_roots = load_pvs_sources(pvs_flow_file)
    netw = u.function_space().mesh()
    coords = netw.coordinates()
    sources = []
    EXIT = 2
    artery_radii = as_P0_function(artery_radii)
    root_idx = np.where(artery_roots.array()>0)
    terminal_points = coords[root_idx,:].squeeze()
    n = FacetNormal(netw)
    ut = TrialFunction(artery_radii.function_space())
    v = TestFunction(artery_radii.function_space())
    un = Function(artery_radii.function_space())
    ds = Measure("ds", netw, subdomain_data=artery_roots)
    solve(v*ut*ds == dot(u, n)*v*ds, un)
    vol_inside = 0
    vol_outside = 0
    from vtk_adapter import create_vtk_structures
    topology, cell_types, x = create_vtk_structures(sas_mesh_func.function_space())
    grid = pv.UnstructuredGrid(topology, cell_types, x)
    grid.point_data["original_ids"] = np.arange(grid.n_points)

    cc_points = grid.cell_centers().points
    tree = KDTree(cc_points)

    _, mapped_ids = tree.query(terminal_points)
    for tp,mi in zip(terminal_points, mapped_ids):
        u_tp = u(*tp)
        u_mag = np.linalg.norm(u_tp)
        assert np.isclose(u_mag, abs(un(*tp)))
        r_art = artery_radii(*tp)
        r_pvs = 2*r_art
        A = np.pi*(r_pvs**2 - r_art**2)
        tip_loc = tp
        try:
            sas_mesh_func(*tp)
        except RuntimeError:
            vol_outside += u_mag * A
            tip_loc = cc_points[mi,:]
            #continue
        vol_inside += u_mag * A
        sources.append(FluxPointSource(x=tip_loc, 
                            normal=u_tp / u_mag,
                            radius=r_pvs,
                            value=un(*tp)))


    pv.PointSet([s.x for s in sources]).cast_to_unstructured_grid().save("tps.vtk")

    print(f"u inside: {vol_inside*m3s_to_mlday:.2f} ml/day")
    print(f"u outside (not applied): {vol_outside*m3s_to_mlday:.2f} ml/day")

    pvs_radii = artery_radii * 2
    A = np.pi*(pvs_radii**2 - artery_radii**2)
    abs_tot_vol = assemble(A*sqrt(dot(u, u))*ds)
    #assert np.isclose(abs_tot_vol*m3s_to_mlday, (vol_outside + vol_inside)*m3s_to_mlday) # fails!
    #assert (artery_roots.array()==1).sum() == assemble(1*ds(1))
    print(abs_tot_vol*m3s_to_mlday)
    
    point_mass_sources = []
    for source in sources:
        point_mass_sources.append((Point(source.x), - source.value* source.radius**2*np.pi))    

    return sources, point_mass_sources

def directsolve(a, L, bcs, W, pointsource=None):
    wh = ii_Function(W)

    A, b = map(ii_assemble, (a, L))
    if pointsource is not None:
        pointsource.apply(b[1])
    A, b = apply_bc(A, b, bcs=bcs)
    A, b = map(monolithic, (A, b))

    solver = PETScLUSolver("mumps")
    solver.solve(A, wh.vector(), b)

    return wh

PETScOptions.set("mat_mumps_icntl_4", 3)  # mumps verbosity
#PETScOptions.set("mat_mumps_icntl_24", 1)  # null pivot detection
#PETScOptions.set("mat_mumps_cntl_7", 1e-8)  # BLR eps
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


def get_normal_func(mesh):
    n = FacetNormal(mesh)
    V  = VectorFunctionSpace(mesh, 'DG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    ds0 = Measure("ds", metadata={'quadrature_degree': 1})
    a = inner(u,v)*ds0
    l = inner(n, v)*ds0
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L, "cg", "jacobi")
    np.testing.assert_allclose(assemble(1*ds(domain=mesh)),
                               assemble(inner(n, nh)*ds), rtol=1e-10, atol=1e-12)
    return nh

def get_inflow_bcs(V, ds, inflow_bcs):
    bcs = []
    n = FacetNormal(V.mesh())
    for bids, infl in inflow_bcs:
        if infl == 0:
            infl_func = Constant((0,0,0))
        else:
            area = 0
            for bid in bids: area += assemble(1*ds(bid))
            infl_func = get_normal_func(V.mesh())
            A = area
            infl_func.vector()[:] *= -eval(infl)

        actual_inflow  = sum(assemble(inner(-n, infl_func)*ds(bid)) for bid in bids)
        infl_num = float(str(infl).translate({ord(i): None for i in 'A*/'}))

        np.testing.assert_allclose(actual_inflow, infl_num, rtol=1e-10, atol=1e-12)
        for bid in bids:
            bcs += [DirichletBC(V, infl_func, 
                                ds.subdomain_data(), bid)]
    return bcs

def get_BDM_problem(mesh, ds, order_k, mu, resistance_bcs, 
                    no_slip_bcs, inflow_bcs, sources):
    
    cell = mesh.ufl_cell()  
    Velm = FiniteElement('BDM', cell, order_k)
    Qelm = FiniteElement('Discontinuous Lagrange',
                          cell, order_k  - 1) 

    V = FunctionSpace(mesh, Velm)
    Q = FunctionSpace(mesh, Qelm)
    W = (V, Q)

    u, p = map(TrialFunction, W) 
    v, q = map(TestFunction, W)

    n = FacetNormal(mesh)
    penalty = Constant(10*order_k)
    CellDiameter = CentroidDistance
    hF = CellDiameter(mesh)
    mu = Constant(float(mu))
    f = Constant([0]*mesh.geometric_dimension())

    a, L = block_form(W, 2), block_form(W, 1)
    a[0][0] = inner(mu*grad(u), grad(v))*dx
    a[0][1] = - inner(p, div(v))*dx
    a[1][0] = - inner(q, div(u))*dx
    
    a[0][0] += Stabilization(mesh, u,v, mu, penalty)
    L[0] = inner(f, v)*dx

    for bid, R in resistance_bcs:
        a[0][0] += inner(Constant(float(R))*dot(u,n), dot(v,n))*ds(bid)

    for bid in no_slip_bcs:
        a[0][0] += -inner(dot(mu*grad(v), n), Tangent(u, n))*ds(bid) \
                - inner(dot(mu*grad(u), n), Tangent(v, n))*ds(bid) \
                + 2*mu*(penalty/hF)*inner(Tangent(u, n), Tangent(v, n))*ds(bid)

    Vbcs = get_inflow_bcs(V, ds, inflow_bcs)
    Qbcs = []
    Wbcs = [Vbcs, Qbcs]
    return a, L, Wbcs, W

def get_BDM_problem_sources(mesh, ds, order_k, mu, resistance_bcs, 
                    no_slip_bcs, inflow_bcs, sources):
    
    cell = mesh.ufl_cell()  
    Velm = FiniteElement('BDM', cell, order_k)
    Qelm = FiniteElement('Discontinuous Lagrange',
                          cell, order_k  - 1) 

    #   Multiplier will like on the auxiliary mesh
    line_facet_f = point_source_mesh(sources, h=mesh.hmin())
    line_mesh = line_facet_f.mesh()
    #from IPython import embed; embed()
    ds_ = Measure('ds', domain=line_mesh, subdomain_data=line_facet_f)
    
    dimM = len(sources)
    M = VectorFunctionSpace(line_mesh, 'R', 0, dimM)

    V = FunctionSpace(mesh, Velm)
    Q = FunctionSpace(mesh, Qelm)
    W = (V, Q, M)

    u, p, l = map(TrialFunction, W) 
    v, q, m = map(TestFunction, W)

    n = FacetNormal(mesh)
    penalty = Constant(10*order_k)
    CellDiameter = CentroidDistance
    hF = CellDiameter(mesh)
    mu = Constant(float(mu))
    f = Constant([0]*mesh.geometric_dimension())

    a, L = block_form(W, 2), block_form(W, 1)
    a[0][0] = inner(mu*grad(u), grad(v))*dx
    a[0][1] = - inner(p, div(v))*dx
    a[1][0] = - inner(q, div(u))*dx
    
    a[0][0] += Stabilization(mesh, u,v, mu, penalty)
    L[0] = inner(f, v)*dx

    for (k, source) in enumerate(sources, 1):
        disk = Disk(radius=source.radius, degree=20, quad_scheme='simple')
        Pi_u, Pi_v = (Average(arg, line_mesh, disk) for arg in (u, v))        

        n_disk = Constant(source.normal)
        area = Constant(pi*source.radius**2)

        assert assemble(Constant(1)*ds_(k)) > 0
        # n_disk is normal from the Stokes side
        a[2][0] += area*inner(dot(Pi_u, n_disk), m[k-1])*ds_(k)
        a[0][2] += area*inner(dot(Pi_v, n_disk), l[k-1])*ds_(k)

        # Here value is u_D.n_D (darcy flux . normal wrt Darcy domain)
        # Check signs based on what will be the input format
        L[2] += -area*inner(Constant(source.value), m[k-1])*ds_(k)    


    for bid, R in resistance_bcs:
        a[0][0] += inner(Constant(float(R))*dot(u,n), dot(v,n))*ds(bid)

    for bid in no_slip_bcs:
        a[0][0] += -inner(dot(mu*grad(v), n), Tangent(u, n))*ds(bid) \
                - inner(dot(mu*grad(u), n), Tangent(v, n))*ds(bid) \
                + 2*mu*(penalty/hF)*inner(Tangent(u, n), Tangent(v, n))*ds(bid)

    Vbcs = get_inflow_bcs(V, ds, inflow_bcs)
    Qbcs = []
    Lbcs = []
    Wbcs = [Vbcs, Qbcs, Lbcs]
    
    return a, L, Wbcs, W


def compute_sas_flow(configfile : str):
    config = read_config(configfile)
    modelname = Path(configfile).stem
    if config.get("out-of-core-factorization", 0):
        PETScOptions.set("mat_mumps_icntl_22", )
    else:
        PETScOptions.set("mat_mumps_icntl_35", 1)  # BLR feature

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


    sources, point_mass_sources = get_pvs_point_sources(Function(FunctionSpace(mesh, "CG", 1)))
    sources = []
    print(f"number of sources {len(sources)}")
    if config["discretization"] == "BDM":
        a, L, bcs, W = get_BDM_problem(mesh, ds, order_k, config["mu"],
                                       config["resistance_bcs"],
                                       config["no_slip_bcs"],
                                       config["inflow_bcs"],
                                       sources)

    else: print("choose BDM discretization!"); exit()

    pointsource = PointSource(W[1], point_mass_sources)

    wh = directsolve(a, L, bcs, W, pointsource)
    uh, ph = wh
    n = FacetNormal(mesh)
    for bids, infl in config["inflow_bcs"]:
        actual_inflow  = sum(assemble(inner(-n, uh)*ds(bid)) for bid in bids)
        print(infl)
        infl_num = float(str(infl).translate({ord(i): None for i in 'A*/'}))
        print(actual_inflow, infl_num)
        np.testing.assert_allclose(actual_inflow, infl_num, rtol=1e-4, atol=1e-12)

    uhdg = interpolate(uh, VectorFunctionSpace(mesh, "DG", order_k))
    absdiv = assemble(abs(div(uhdg))*dx)
    print(f"absdiv = {absdiv*m3s_to_mlday}")
    if ph.ufl_element().family() == 'Discontinuous Lagrange' and pointsource is None:
        assert np.isclose(absdiv, 0)
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
    print(f"tot_outflow = {assemble(dot(uh, n)*ds)*m3s_to_mlday}")
    print(f"LV_outflow  = {assemble(dot(uh, n)*ds(LV_INTERF_ID))*m3s_to_mlday}")
    print(f"UPPER_SKULL_outflow = {assemble(dot(uh, n)*ds(UPPER_SKULL_ID))*m3s_to_mlday}")
    print(f"abs PVS flow: {np.sum([abs(s.value) * s.radius**2 * np.pi for s in sources])*m3s_to_mlday}")
    print(f"net PVS inflow: {np.sum([s.value * s.radius**2 * np.pi for s in sources])*m3s_to_mlday}")

    
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
