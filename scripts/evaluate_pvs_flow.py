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
from plotting_utils import minmax
import ufl_legacy as ufl
m2mm = 1e3
m2mum = 1e6
Pa2mPa = 1e3

def sig_to_space(sig, mesh):
    elem = eval(sig.replace("Cell", "ufl.Cell"))
    return FunctionSpace(mesh, elem)

def compute_pvs_flow(pvs_flow_file):

    mesh = Mesh()

    with HDF5File(MPI.comm_world, pvs_flow_file,'r') as f:
        f.read(mesh, "mesh", False)
        p_sig = f.attributes("/p").to_dict()["signature"]
        umag_sig = f.attributes("/umag").to_dict()["signature"]
        r_sig = f.attributes("/radii").to_dict()["signature"]
        p_space, umag_space, r_space = map(lambda sig: sig_to_space(sig, mesh), [p_sig, umag_sig, r_sig])
        p, uh_mag, artery_radii = map(Function, [p_space, umag_space, r_space])
        for u,n in zip([p, uh_mag, artery_radii], ["p", "umag", "radii"]):
            f.read(u, n)
    modelstr = pvs_flow_file.split('/')[1:-1]
    model = modelstr[-1]
    plot_dir = f"plots/{'/'.join(modelstr)}"
    os.makedirs(plot_dir, exist_ok=True)
    segments, segids ,_ = color_branches(mesh)
    dxs = Measure("dx", mesh, subdomain_data=segments)

    seglengths = np.array([assemble(1*dxs(si)) for si in segids])
    assert (seglengths > 0).all()
    umag = np.array([assemble(m2mum*uh_mag*dxs(si))/ l for si,l in zip(segids, seglengths)])
    umagabs = abs(umag)


    set_plotting_defaults()
    celllengths = project(CellDiameter(mesh), FunctionSpace(mesh, "DG", 0)).vector()[:]

    for variant, uvals, lengths in zip(["cell", "seg"],
                                       [uh_mag.vector().get_local()*m2mum, umag],
                                       [celllengths , seglengths]):
        uavg = np.average(uvals, weights=lengths)
        umax = abs(uvals).max()
        umed = np.median(abs(uvals))
        range = np.round(minmax([uvals], percentile=99.9), 1)
        nbins = 50
        fig, ax = plt.subplots(figsize=(4,3))
        counts, bins, containers = plt.hist(uvals, density=False, bins=nbins, histtype="bar",
                range=range, edgecolor='black', linewidth=0.4, 
                        weights=lengths / lengths.sum())
        pos_share = max([0, range[1] / (range[1] - range[0])])
        print(pos_share)
        for bar, ca in zip(containers, sns.color_palette("Purples_r",
                                                         n_colors=int((1 - pos_share) * nbins)) + 
                                       sns.color_palette("Greens",
                                                       n_colors=int(pos_share * nbins))):
            bar.set_facecolor(ca)
        plt.axvline(uavg, color="black", linestyle=":", label=f"avg: {uavg:.2f} µm/s")
        #plt.axvline(umed, color="black", linestyle="-.", label=f"median: {umed:.2f} µm/s")
        plt.axvline(umax, color="black", linestyle="--", label=f"max: {umax:.2f} µm/s")
        plt.xlabel("PVS flow velocity (µm/s)")
        plt.ylabel("frequency")
        plt.tight_layout()
        plt.xlim(range)
        plt.legend(frameon=False)
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        plt.savefig(f"{plot_dir}/{model}_velocity_histo_{variant}.png", dpi=300)

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
        pressures.append(assemble(Pa2mPa*p*dxs(si)) / length)
        velocities.append(assemble(m2mum*uh_mag*dxs(si))/ length)
        radii.append(assemble(m2mm*artery_radii*dxs(si))/ length)
        lengths.append(length*m2mm)
    df = pd.DataFrame({"p":pressures, "u":velocities, "loc":artlabels, "L":lengths,
                        "R":radii})

    plt.figure(figsize=(4,3))
    ax = sns.barplot(df, x="loc", y="p",)# palette="crest")
    ax.set_xlabel("")
    ax.set_ylabel("pressure (mPa)")
    plt.xticks(rotation=45,ha='right', rotation_mode='anchor', size="10");plt.tight_layout()
    plt.savefig(f"{plot_dir}/{model}_pressure.png")

    plt.figure(figsize=(4,3))
    ax = sns.barplot(df, x="loc", y="u",)# palette="flare")
    ax.set_xlabel("")
    ax.set_ylabel("flow velocity (μm/s)")
    plt.xticks(rotation=45,ha='right', rotation_mode='anchor', size="10");plt.tight_layout()
    plt.savefig(f"{plot_dir}/{model}_velocity.png")

    plt.figure(figsize=(4,3))
    ax = sns.barplot(df, x="loc", y="R",)# palette="blend:#7AB,#EDA")
    ax.set_xlabel("")
    ax.set_ylabel("radius (mm)")
    plt.xticks(rotation=45,ha='right', rotation_mode='anchor', size="10");plt.tight_layout()
    plt.savefig(f"{plot_dir}/{model}_radius.png")

    plt.figure(figsize=(4,3))
    ax = sns.barplot(df, x="loc", y="L",)# palette="blend:#7AB,#EDA")
    ax.set_xlabel("")
    ax.set_ylabel("segment length (mm)")
    plt.xticks(rotation=45,ha='right', rotation_mode='anchor', size="10");plt.tight_layout()
    plt.savefig(f"{plot_dir}/{model}_length.png")

    length = assemble(1*dx(domain=mesh))
    umean = assemble(uh_mag*dx) / length
    umax = norm(uh_mag.vector(), "linf")
    metrics = dict(pvsumean=umean,pvsumax=umax)
    with open(f'{plot_dir}/pvs_metrics.yml', 'w') as outfile:
        yaml.dump(metrics, outfile, default_flow_style=False)

if __name__ == "__main__":
    typer.run(compute_pvs_flow)
