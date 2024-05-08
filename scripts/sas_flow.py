from dolfin import *
from solver import * 
import os
import numpy as np
from xii import *
from petsc4py import PETSc
import sympy as sp
import typer

PETScOptions.set("mat_mumps_icntl_4", 3)  # mumps verbosity
#PETScOptions.set("mat_mumps_icntl_24", 1)  # null pivot detection
PETScOptions.set("mat_mumps_icntl_35", 1)  # BLR feature
#PETScOptions.set("mat_mumps_icntl_32", 1)  # forward elimination during solve (useful, but not passed on by petsc)
PETScOptions.set("mat_mumps_icntl_22", 1)  # out-of-core to reduce memory
#PETScOptions.set("mat_mumps_icntl_11", 1)  # error analysis
#PETScOptions.set("mat_mumps_icntl_25", 2)  # turn on null space basis


def cell_to_facet_meshfunc(cellfunc):
    facetfct = MeshFunction('size_t', cellfunc.mesh(), cellfunc.dim() - 1)
    facetfct.set_all(0)
    #mesh.init(dim - 1, dim) # relates facets to cells
    for f in facets(cellfunc.mesh()):
        if 1 in cellfunc.array()[f.entities(cellfunc.dim())]: # one of the neighbouring cells of the facet is part of region 1
            facetfct[f.index()] = 1
    return facetfct


def compute_sas_flow(meshname : str):
# get mesh 
    sas = Mesh()
    with XDMFFile(f'mesh/{meshname}/volmesh/colors.xdmf') as f:
        f.read(sas)
        gdim = sas.geometric_dimension()
        sas_components = MeshFunction('size_t', sas, gdim, 0)
        f.read(sas_components, 'sas_components')
        sas.scale(1e-3)  # scale from mm to m

    sas_outer = EmbeddedMesh(sas_components, [1,3,4]) 

    # create boundary markers 
    boundary_markers = MeshFunction("size_t", sas, sas.topology().dim() - 1)
    boundary_markers.set_all(0)
    # Sub domain for efflux route (mark whole boundary of the full domain) 
    class Efflux(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[2] > 0.1
    efflux = Efflux()
    efflux.mark(boundary_markers, 1)
    
    # translate markers to the sas outer mesh 
    boundary_markers_outer = MeshFunction("size_t", sas_outer, sas_outer.topology().dim() - 1, 0)
    DomainBoundary().mark(boundary_markers_outer, 2)
    sas_outer.translate_markers(boundary_markers, (1, ), marker_f=boundary_markers_outer)

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
    R = Constant(1e-2) # 1e-5 Pa/(mm s)
    cp1_midpoint = [0.128, 0.229, 0.192] # found in paraview
    cp2_midpoint = [0.2, 0.229, 0.192] # found in paraview

    #g1 =  Expression("exp( - ((x[0] - m0)*(x[0] - m0) + (x[1] - m1)*(x[1] - m1) + (x[2] - m2)*(x[2] - m2)) / (sigma*sigma))",
    #                    m0=cp1_midpoint[0], m1=cp1_midpoint[1],
    #                    m2=cp1_midpoint[2], sigma=0.01, degree=3)
    #g2 =  Expression("exp( - ((x[0] - m0)*(x[0] - m0) + (x[1] - m1)*(x[1] - m1) + (x[2] - m2)*(x[2] - m2)) / (sigma*sigma))",
    #                    m0=cp2_midpoint[0], m1=cp2_midpoint[1],
    #                    m2=cp2_midpoint[2], sigma=0.01, degree=3)

    LV_marker = as_P0_function(sas_components)
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

    uh, ph = wh.split(deepcopy=True)[:2]
    
    CG2 = VectorFunctionSpace(sas_outer, "CG", 2)
    uh2 = interpolate(uh, CG2)
    with XDMFFile(f'results/csf_flow/{meshname}/csf_vis.xdmf') as xdmf:
        xdmf.write_checkpoint(uh2, "velocity")
        xdmf.write_checkpoint(ph, "pressure", append=True)

    # extend the solution by zero on the whole domain:
    CG3 = VectorFunctionSpace(sas, "CG", 3)
    uh.set_allow_extrapolation(True)
    uh_global = interpolate(uh, CG3)
    sas_fct_func = cell_to_facet_meshfunc(sas_components)
    bc = DirichletBC(CG3, Constant([0]*gdim), sas_fct_func, 0)
    bc.apply(uh_global.vector())

    with XDMFFile(f'results/csf_flow/{meshname}/csf_v.xdmf') as xdmf:
        xdmf.write_checkpoint(uh_global, "velocity")

    with XDMFFile(f'results/csf_flow/{meshname}/csf_p.xdmf') as xdmf:
        xdmf.write_checkpoint(ph, "pressure")

if __name__ == "__main__":
    typer.run(compute_sas_flow)