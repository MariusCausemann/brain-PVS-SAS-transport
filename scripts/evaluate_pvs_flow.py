from dolfin import *
import numpy as np
import os
from test_map_on_global_coords_shift import map_kdtree, map_dg_on_global
from subdomain_ids import CSFID, V34ID, PARID, LVID, CSFNOFLOWID
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
from cmap import Colormap 
from vtk_adapter import create_vtk_structures
from extract_vessels import get_tubes
import ufl_legacy as ufl

m2mm = 1e3
m2mum = 1e6
Pa2mPa = 1e3

def plot_pvs_velocity(mesh, uh_mag, artery_radii, filename):
    CG1 = FunctionSpace(mesh, "CG", 1)
    bar_args=dict(title="velocity (μm/s)", vertical=False,
                    height=0.07, width=0.6, position_x=0.2,
                    position_y=0.0, title_font_size=28, font_family="arial",
                    label_font_size=28, fmt="%.0f")
    topology, cell_types, x = create_vtk_structures(CG1)
    grid = pv.UnstructuredGrid(topology, cell_types, x)
    grid["u"] = uh_mag.vector()[:] *m2mum
    umax = abs(grid["u"]).max() 
    grid["radius"] = artery_radii.vector()[:]
    tubes = get_tubes(grid.ctp())
    pl = pv.Plotter()
    pl = pv.Plotter(off_screen=True, window_size=(1200, 800))
    pl.add_mesh(tubes,scalars="u", cmap="curl", clim=(-umax, umax),
                scalar_bar_args=bar_args)
    pl.camera_position = 'xz'
    #pl.camera.roll += 5
    pl.camera.azimuth += 165
    pl.camera.elevation += 10
    pl.camera.zoom(1.85)
    pl.camera.focal_point = np.array(pl.camera.focal_point) + np.array([0,0,-0.008])
    img = pl.screenshot(filename, transparent_background=True)

def sig_to_space(sig, mesh):
    elem = eval(sig.replace("Cell", "ufl.Cell"))
    return FunctionSpace(mesh, elem)

def compute_pvs_flow(pvs_flow_file):
    artcolor = "#94221F"
    color = "#95AAD3"
    mesh = Mesh()

    with HDF5File(MPI.comm_world, pvs_flow_file,'r') as f:
        f.read(mesh, "mesh", False)
        p_sig = f.attributes("/p").to_dict()["signature"]
        umag_sig = f.attributes("/umag").to_dict()["signature"]
        r_sig = f.attributes("/radii").to_dict()["signature"]
        p_space, umag_space, r_space = map(lambda sig: sig_to_space(sig, mesh), [p_sig, umag_sig, r_sig])
        p, uh_mag, artery_radii, pvs_radii = map(Function, [p_space, umag_space, r_space, r_space])
        for u,n in zip([p, uh_mag, artery_radii, pvs_radii], 
                       ["p", "umag", "radii", "pvs_radii"]):
            f.read(u, n)
    modelstr = pvs_flow_file.split('/')[1:-1]
    model = modelstr[-1]
    plot_dir = f"plots/{'/'.join(modelstr)}"
    os.makedirs(plot_dir, exist_ok=True)

    plot_pvs_velocity(mesh, uh_mag, pvs_radii,
                      plot_dir + "/vel3D.png")

    segments, segids ,_ = color_branches(mesh)
    dxs = Measure("dx", mesh, subdomain_data=segments)
    segs = segments.array().copy()
    seglengths, umag = [], []
    for si in segids:
        segments.array()[:] = np.where(segs==si, 1, 0)
        length = assemble(1*dxs(1))
        mag = assemble(m2mum*uh_mag*dxs(1))/ length
        umag.append(mag)
        seglengths.append(length)
    segments.array()[:] = segs
    seglengths = np.array(seglengths)
    umag = np.array(umag)
    assert (seglengths > 0).all()
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
        fig, ax = plt.subplots(figsize=(5,2.8))
        counts, bins, containers = plt.hist(uvals, density=False, bins=nbins, histtype="bar",
                range=range, edgecolor='black', linewidth=0.4, 
                        weights=lengths / lengths.sum())
        pos_share = max([0, range[1] / (range[1] - range[0])])
        print(pos_share)
        for bar, ca in zip(containers, list(Colormap("balance_blue").iter_colors(int((1 - pos_share) * nbins))) + 
                                       list(Colormap("curl_pink_r").iter_colors(int(pos_share * nbins)))):
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
        plt.savefig(f"{plot_dir}/{model}_velocity_histo_{variant}.png", dpi=300,
                    transparent=True)

    points = np.array([c for n,c in pointlabels])
    cellidx = map_kdtree(FunctionSpace(mesh, "DG", 0).tabulate_dof_coordinates(),
                        points)

    artlabels = [n for n,c in pointlabels]
    artsegids = [int(segments.array()[i]) for i in cellidx]

    # plot pressure, velocity and radius at main arteries
    pressures, velocities, radii, lengths = [], [], [], []
    segs = segments.array().copy()
    for si in artsegids:
        segments.array()[:] = np.where(segs==si, 1, 0)
        length = assemble(1*dxs(1))
        assert length > 0
        pressures.append(assemble(Pa2mPa*p*dxs(1)) / length)
        velocities.append(assemble(m2mum*uh_mag*dxs(1))/ length)
        radii.append(assemble(m2mm*artery_radii*dxs(1))/ length)
        lengths.append(length*m2mm)
    df = pd.DataFrame({"p":pressures, "u":velocities, "loc":artlabels, "L":lengths,
                        "R":radii})

    plt.figure(figsize=(4,3))
    ax = sns.barplot(df, x="loc", y="p", color=color)# palette="crest")
    ax.set_xlabel("")
    ax.set_ylabel("pressure (mPa)")
    plt.xticks(rotation=45,ha='right', rotation_mode='anchor', size="10");plt.tight_layout()
    plt.subplots_adjust(left=0.2, bottom=0.3, right=0.98)
    plt.savefig(f"{plot_dir}/{model}_pressure.png", transparent=True)

    plt.figure(figsize=(4,3))
    ax = sns.barplot(df, x="loc", y="u",color=color)# palette="flare")
    ax.set_xlabel("")
    ax.set_ylabel("flow velocity (μm/s)")
    plt.xticks(rotation=45,ha='right', rotation_mode='anchor', size="10");plt.tight_layout()
    plt.subplots_adjust(left=0.2, bottom=0.3, right=0.98)
    plt.savefig(f"{plot_dir}/{model}_velocity.png", transparent=True)

    plt.figure(figsize=(4,3))
    ax = sns.barplot(df, x="loc", y="R",color=artcolor)# palette="blend:#7AB,#EDA")
    ax.set_xlabel("")
    ax.set_ylabel("radius (mm)")
    plt.xticks(rotation=45,ha='right', rotation_mode='anchor', size="10");plt.tight_layout()
    plt.subplots_adjust(left=0.2, bottom=0.3, right=0.98)
    plt.savefig(f"{plot_dir}/{model}_radius.png", transparent=True)

    plt.figure(figsize=(4,3))
    ax = sns.barplot(df, x="loc", y="L",color=color)# palette="blend:#7AB,#EDA")
    ax.set_xlabel("")
    ax.set_ylabel("segment length (mm)")
    plt.xticks(rotation=45,ha='right', rotation_mode='anchor', size="10");plt.tight_layout()
    plt.subplots_adjust(left=0.2, bottom=0.3, right=0.98)
    plt.savefig(f"{plot_dir}/{model}_length.png", transparent=True)

    length = assemble(1*dx(domain=mesh))
    umean = assemble(uh_mag*dx) / length
    umax = norm(uh_mag.vector(), "linf")
    metrics = dict(pvsumean=umean,pvsumax=umax)
    with open(f'{plot_dir}/pvs_metrics.yml', 'w') as outfile:
        yaml.dump(metrics, outfile, default_flow_style=False)

if __name__ == "__main__":
    typer.run(compute_pvs_flow)
