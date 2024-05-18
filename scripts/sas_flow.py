from dolfin import *
from solver import * 
import os
import numpy as np
from xii import *
from petsc4py import PETSc
import sympy as sp
import typer
import numpy_indexed as npi
from IPython import embed

PETScOptions.set("mat_mumps_icntl_4", 3)  # mumps verbosity
#PETScOptions.set("mat_mumps_icntl_24", 1)  # null pivot detection
PETScOptions.set("mat_mumps_icntl_35", 1)  # BLR feature
#PETScOptions.set("mat_mumps_icntl_32", 1)  # forward elimination during solve (useful, but not passed on by petsc)
PETScOptions.set("mat_mumps_icntl_22", 1)  # out-of-core to reduce memory
#PETScOptions.set("mat_mumps_icntl_11", 1)  # error analysis
#PETScOptions.set("mat_mumps_icntl_25", 2)  # turn on null space basis

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

def cell_to_facet_meshfunc(cellfunc, label):
    facetfct = MeshFunction('size_t', cellfunc.mesh(), cellfunc.dim() - 1)
    facetfct.set_all(0)
    #mesh.init(dim - 1, dim) # relates facets to cells
    for f in facets(cellfunc.mesh()):
        if label in cellfunc.array()[f.entities(cellfunc.dim())]: # one of the neighbouring cells of the facet is part of region 1
            facetfct[f.index()] = 1
    return facetfct

def map_on_global(uh, uh_global):
    # map the solution back on the whole domain
    cell_map = uh.function_space().mesh().parent_entity_map[0][3]
    for child, parent in cell_map.items():
        child_dofs = uh.function_space().dofmap().cell_dofs(child)
        parent_dofs = uh_global.function_space().dofmap().cell_dofs(parent)
        uh_global.vector().vec().array_w[parent_dofs] = uh.vector().vec().array_r[child_dofs]


def interpolate_on_global(uh, uh_global, label):

    uh.set_allow_extrapolation(True)
    uh_global.interpolate(uh)
    dim = uh.geometric_dimension()
    lblfct_par = cell_to_facet_meshfunc(label, PARID)
    lblfct_noflow = cell_to_facet_meshfunc(label, CSFNOFLOWID)
    bcmagglobpar = DirichletBC(uh_global.function_space(), Constant([0]*dim), lblfct_par, 1)
    bcmagglobnoflow = DirichletBC(uh_global.function_space(), Constant([0]*dim), lblfct_noflow, 1)
    bcmagglobpar.apply(uh_global.vector())
    bcmagglobnoflow.apply(uh_global.vector())

def compute_sas_flow(meshname : str):
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
    class Efflux(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary
    efflux = Efflux()
    efflux.mark(boundary_markers, 1)
    
    # translate markers to the sas outer mesh 
    boundary_markers_outer = MeshFunction("size_t", sas_outer, sas_outer.topology().dim() - 1, 0)
    DomainBoundary().mark(boundary_markers_outer, 2)
    #sas_outer.translate_markers(boundary_markers, (1, ), marker_f=boundary_markers_outer)

    File("fpp.pvd") << boundary_markers_outer

    ds = Measure("ds", domain=sas_outer, subdomain_data=boundary_markers_outer)

    # Define function spaces for velocity and pressure
    cell = sas_outer.ufl_cell()  
    Velm = VectorElement('Lagrange', cell, 3)
    Qelm = FiniteElement('Lagrange', cell, 2) 
    Q = FunctionSpace(sas_outer, Qelm)
    W = FunctionSpace(sas_outer, Velm * Qelm)

    u , p = TrialFunctions(W) 
    v , q = TestFunctions(W)
    n = FacetNormal(sas_outer)

    mu = Constant(0.7e-3) # units need to be checked 
    R = Constant(1e4) # 1e-5 Pa/(mm s)

    LV_marker = as_P0_function(label)
    LV_marker.vector()[:] = LV_marker.vector()[:] == 3
    LV_marker_outer = interpolate(LV_marker, FunctionSpace(sas_outer, "DG", 0))
    g1 = g2 = LV_marker_outer

    gtot = assemble((g1 + g2 )*dx(domain=sas_outer))

    total_production = 0.63e-3/(60*60*24) # 0.63 L / day
    g = (g1+ g2)* Constant(total_production / gtot)

    f = Constant([0]*gdim)
    g_L_per_day = 1e3 *(60*60*24) * assemble(g*dx(domain=sas_outer))
    print(f"production rate: {g_L_per_day} L/day")
    assert np.isclose(g_L_per_day, 0.63)
    gfunc = project(g, Q)
    File("g.pvd") << gfunc
    a = (inner(2*mu*sym(grad(u)), sym(grad(v)))*dx - inner(p, div(v))*dx
            -inner(q, div(u))*dx + inner(R*dot(u,n), dot(v,n))*ds(1)) 
    L = inner(f, v)*dx - g*q*dx
    bcs = [DirichletBC(W.sub(0), Constant((0, 0, 0)), ds.subdomain_data(), 2)] 

    wh = Function(W)

    solve(a==L, wh, bcs=bcs, solver_parameters={"linear_solver":"mumps"})

    uh, ph = wh.split(deepcopy=True)[:]
    
    # project to CG2 and write for visualization
    CG2 = VectorFunctionSpace(sas_outer, "CG", 2)
    uh2 = project(uh, CG2, solver_type="cg", preconditioner_type="hypre_amg")
    with XDMFFile(f'results/csf_flow/{meshname}/csf_vis_v_{R}.xdmf') as xdmf:
        xdmf.write_checkpoint(uh2, "velocity")
    with XDMFFile(f'results/csf_flow/{meshname}/csf_vis_p_{R}.xdmf') as xdmf:
        xdmf.write_checkpoint(ph, "pressure")

    CG3 = VectorFunctionSpace(sas, "CG", 3)
    uh_global = Function(CG3)
    interpolate_on_global(uh, uh_global, label)
    dxglob = Measure("dx", sas, subdomain_data=label)

    assert np.isclose(assemble(inner(uh_global, uh_global)*dxglob(PARID)), 0)
    assert np.isclose(assemble(inner(uh_global, uh_global)*dxglob(CSFNOFLOWID)), 0)
    assert np.isclose(assemble(inner(uh, uh)*dx), assemble(inner(uh_global, uh_global)*dxglob))

    divu_global = assemble(div(uh_global)*div(uh_global)*dxglob)
    divu = assemble(div(uh)*div(uh)*dx)

    assert np.isclose(divu, divu_global)

    with XDMFFile(f'results/csf_flow/{meshname}/csf_v.xdmf') as xdmf:
        xdmf.write_checkpoint(uh_global, "velocity")

    with XDMFFile(f'results/csf_flow/{meshname}/csf_p.xdmf') as xdmf:
        xdmf.write_checkpoint(ph, "pressure")

    umag = project(sqrt(inner(uh, uh)), FunctionSpace(sas_outer, "CG", 3),
                   solver_type="cg", preconditioner_type="hypre_amg")
    umean = assemble(umag*dx) / assemble(1*dx(domain=sas_outer))

    print(f"div u = {divu}")
    print(f"u max = {umag.vector().max()}")
    print(f"u umean = {umean}")

if __name__ == "__main__":
    typer.run(compute_sas_flow)