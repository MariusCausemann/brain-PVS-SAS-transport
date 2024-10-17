from dolfin import *
import xii
from solver import read_vtk_network, as_P0_function, get_mesh
import numpy as np
import os
from pvs_flow import pvs_flow_system
from test_map_on_global_coords_shift import map_kdtree, map_dg_on_global
from generate_synthseg_mesh import CSFID, V34ID, PARID, LVID, CSFNOFLOWID
import typer
from plotting_utils import read_config, set_plotting_defaults
import seaborn as sns
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from label_arteries import pointlabels
import pandas as pd
import yaml
from branch_marking import color_branches

m2mm = 1e3
m2mum = 1e6
Pa2mPa = 1e3

def map_pressure_on_network(p, pnet):
    idx = map_kdtree(p.function_space().tabulate_dof_coordinates(), 
                     pnet.function_space().tabulate_dof_coordinates())
    pnet.vector()[:] = p.vector()[idx]


def get_blood_flow_orientation(tau, radius_f, supply_nodes):
    mesh = tau.function_space().mesh()
    a, L, W = pvs_flow_system(radius_f, tau, 2, f=Constant(0))

    supply_marker = MeshFunction("size_t", mesh, 0, 0)
    supply_marker.array()[supply_nodes] = 1

    bcout = DirichletBC(W.sub(1), Constant(0), "on_boundary")
    bcin = DirichletBC(W.sub(1), Constant(1), supply_marker, 1)

    A, b = assemble_system(a, L, [bcout, bcin])
    wh = Function(W)
    solve(A, wh.vector(), b)
    uh_mag, ph = wh.split(deepcopy=True)
    uh_sign = Function(uh_mag.function_space())
    uh_sign.vector()[:] = np.where(uh_mag.vector()[:] >=0, 1, -1)
    return uh_sign
# --------------------------------------------------------------------

def compute_pvs_flow(csf_flow_model:  str):

    radius_ratio = 2
    mesh, artery_radii, artery_roots = read_vtk_network("mesh/networks/arteries_smooth.vtk", rescale_mm2m=False)
    pressure_file = f"results/csf_flow/{csf_flow_model}/flow.hdf"
    results_dir = f"results/pvs_flow_prod/{csf_flow_model}/"
    csf_flow_config = read_config(f"configfiles/{csf_flow_model}.yml")
    radius_f = as_P0_function(artery_radii)
    
    pmesh = Mesh()
    with HDF5File(MPI.comm_world, pressure_file,'r') as f:
        f.read(pmesh, "mesh", False)
        p_elem = eval(f.attributes("/pressure").to_dict()["signature"])
        DG = FunctionSpace(pmesh, p_elem)
        p = Function(DG)
        f.read(p, "pressure")

    _, sm = get_mesh(csf_flow_config["mesh"].split("_outer")[0] + ".xdmf")
    sas = xii.EmbeddedMesh(sm, [CSFID]) 
    par = xii.EmbeddedMesh(sm, [PARID, LVID, V34ID, CSFNOFLOWID]) 
    p = interpolate(p, FunctionSpace(sas, "DG", 1))

    p.set_allow_extrapolation(True)
    V = FunctionSpace(par, "CG", 1)
    u, v = TrialFunction(V), TestFunction(V)
    sol = Function(V)
    p_dirichletbc = DirichletBC(V, p, "on_boundary")
    solve(inner(grad(u), grad(v))*dx == Constant(0)*v*dx, sol, bcs=p_dirichletbc,
        solver_parameters={"linear_solver":"cg", "preconditioner":"amg"})

    p_combined = map_dg_on_global(p, parentmesh=sm.mesh())
    p_par = map_dg_on_global(interpolate(sol, FunctionSpace(par, "DG", 1)),
                             parentmesh=sm.mesh())

    p_combined.vector()[:] += p_par.vector()[:]

    p_combined.rename("p","p")


    # Grab the tangent of xii; bottom line is that this is vector valued
    # function on the network describing tangent to the edge. Orientation
    # is arbitrary as long same tau is used throught the code
    tau = xii.TangentCurve(mesh)
    supply_nodes = [8065, 7173, 4085]
    downstream = get_blood_flow_orientation(tau, radius_f, supply_nodes)

    a, L, W = pvs_flow_system(radius_f, downstream*tau, radius_ratio, f=Constant(0))

    pbc = Function(W.sub(1).collapse())
    map_pressure_on_network(p_combined, pbc)
    pbc.rename("p","p")

    bc = DirichletBC(W.sub(1), pbc, "on_boundary")

    A, b = assemble_system(a, L, [bc])

    wh = Function(W)
    solve(A, wh.vector(), b)

    uh_mag, ph = wh.split(deepcopy=True)

    pvs_flow_vec = project(uh_mag*tau*downstream, 
        VectorFunctionSpace(mesh, "DG", 1))
    
    ph.rename("p","p")
    pvs_flow_vec.rename("uh", "uh")
    os.makedirs(results_dir, exist_ok=True)
    with XDMFFile(f'{results_dir}/pvs_flow.xdmf') as xdmf:
        xdmf.write_checkpoint(pvs_flow_vec, "velocity")
    with XDMFFile(f'{results_dir}/pvs_flow_vis.xdmf') as xdmf:
        xdmf.write(pvs_flow_vec, t=0)
        xdmf.write(ph, t=0)
    

    segments, segids ,_ = color_branches(mesh)
    dxs = Measure("dx", mesh, subdomain_data=segments)


    seglengths = np.array([assemble(1*dxs(si)) for si in segids])
    assert (seglengths > 0).all()
    umag = np.array([assemble(m2mum*uh_mag*dxs(si))/ l for si,l in zip(segids, seglengths)])
    umagabs = abs(umag)
    uavg = np.average(umagabs, weights=seglengths)
    umax = umagabs.max()
    umed = np.median(umagabs)

    set_plotting_defaults()

    fig, ax = plt.subplots()
    counts, bins, containers = plt.hist(umag, density=False, bins=50, histtype="bar",
            range=(-0.4, 0.6), edgecolor='black', linewidth=0.4, 
                    weights=seglengths / seglengths.sum())
    for bar, ca in zip(containers, sns.color_palette("Purples_r", n_colors=20) + 
                                   sns.color_palette("Greens", n_colors=30)):
        bar.set_facecolor(ca)
    plt.axvline(uavg, color="black", linestyle=":", label=f"avg: {uavg:.2f} µm/s")
    plt.axvline(umed, color="black", linestyle="-.", label=f"median: {umed:.2f} µm/s")
    plt.axvline(umax, color="black", linestyle="--", label=f"max: {umax:.2f} µm/s")
    plt.xlabel("PVS flow velocity (µm/s)")
    plt.ylabel("frequency")
    plt.tight_layout()
    plt.xlim((-0.4, 0.6))
    plt.legend(frameon=False)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    plt.savefig(f"{results_dir}/velocity_histo.png")

    points = np.array([c for n,c in pointlabels])
    cellidx = map_kdtree(FunctionSpace(mesh, "DG", 0).tabulate_dof_coordinates(),
                        points)


    artlabels = [n for n,c in pointlabels]

    artsegids = [int(segments.array()[i]) for i in cellidx]

    # plot pressure, velocity and radius at main arteries
    pressures, velocities, radii, lengths = [], [], [], []
    for si in artsegids:
        length = assemble(1*dxs(si))
        assert length > 0
        pressures.append(assemble(Pa2mPa*ph*dxs(si)) / length)
        velocities.append(assemble(m2mum*uh_mag*dxs(si))/ length)
        radii.append(assemble(m2mm*as_P0_function(artery_radii)*dxs(si))/ length)
        lengths.append(length*m2mm)
    df = pd.DataFrame({"p":pressures, "u":velocities, "loc":artlabels, "L":lengths,
                        "R":radii}).sort_values(by="loc")

    plt.figure()
    ax = sns.barplot(df, x="loc", y="p", palette="crest")
    ax.set_xlabel("")
    ax.set_ylabel("pressure (mPa)")
    plt.savefig(f"{results_dir}/arteries_labels_pressure.png")

    plt.figure()
    ax = sns.barplot(df, x="loc", y="u", palette="flare")
    ax.set_xlabel("")
    ax.set_ylabel("flow velocity (μm/s)")
    plt.savefig(f"{results_dir}/arteries_labels_velocity.png")

    plt.figure()
    ax = sns.barplot(df, x="loc", y="R", palette="blend:#7AB,#EDA")
    ax.set_xlabel("")
    ax.set_ylabel("radius (mm)")
    plt.savefig(f"{results_dir}/arteries_labels_radius.png")

    plt.figure()
    ax = sns.barplot(df, x="loc", y="L", palette="blend:#7AB,#EDA")
    ax.set_xlabel("")
    ax.set_ylabel("segment length (mm)")
    plt.savefig(f"{results_dir}/arteries_labels_length.png")

    length = assemble(1*dx(domain=mesh))
    umean = assemble(uh_mag*dx) / length
    umax = norm(uh_mag.vector(), "linf")
    metrics = dict(pvsumean=umean,pvsumax=umax)
    with open(f'{results_dir}/pvs_metrics.yml', 'w') as outfile:
        yaml.dump(metrics, outfile, default_flow_style=False)

if __name__ == "__main__":
    typer.run(compute_pvs_flow)
