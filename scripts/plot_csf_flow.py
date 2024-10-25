
import pyvista as pv
import typer
from vtk_adapter import create_vtk_structures
import numpy as np
from fenics import *
import matplotlib 
import k3d.colormaps.paraview_color_maps as pcm
from compute_dispersion_field import alpha
from plotting_utils import set_plotting_defaults, minmax
from generate_synthseg_mesh import CSFID, LVID, V34ID
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import seaborn as sns

m2mum = 1e6

def velocity_histo(sm, u, filename):
    set_plotting_defaults()
    DG0 = sm.function_space()
    umag = project(sqrt(inner(u,u)), DG0,
                   solver_type="cg", preconditioner_type="hypre_amg")
    cellvol = assemble(TestFunction(DG0)*dx)
    ids = [CSFID, LVID, V34ID]
    indices = [sm.vector()[:]==i for i in ids]  
    fig, axes = plt.subplots(ncols=3, figsize=(7,3))
    nbins = 20
    colors = sns.color_palette("crest", n_colors=3)
    for ax, name, idx, col in zip(axes, ["SAS", "LV", "V3 & V4"], indices, colors):
        uvals = umag.vector()[idx] * m2mum
        xlims = (0, np.percentile(uvals, 92))
        ax.hist(uvals, bins=nbins, weights=cellvol[idx] / cellvol[idx].sum(),
                 range=xlims,
                  color=col, alpha=0.7, label=name)
        mean = np.average(uvals, weights=cellvol[idx])
        ax.axvline(mean, 
                   color="black",ymax=0.6, label=f"{mean:.1f} µm/s \n (mean)")
        #ax.axvline(np.max(uvals), ls="dotted",
        #           color="black",ymax=0.6, label=f"{max*1e-3:.1f} mm/s")
        ax.set_xlabel("velocity (µm/s)")
        ax.set_ylabel("frequency")
        plt.tight_layout()
        ax.set_xlim(xlims)
        ax.legend(frameon=False, loc="lower left", alignment="right",
                  borderaxespad=0, handlelength=1.3, bbox_to_anchor=[0.2, 0.7])
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    plt.savefig(filename)

def from_k3d(colorlist):
    cm = np.array(colorlist).reshape(-1, 4)[:, 1:]
    return matplotlib.colors.LinearSegmentedColormap.from_list("", cm)

cardiac_config = dict(pcmap="balance", 
                    vcmap="turbo")
prod_config = dict(pcmap=from_k3d([1]*4 + pcm.Blue___Green___Orange), 
                    vcmap="inferno")

def plot_csf_flow(dirname: str):
    filename = f"{dirname}/flow.hdf"

    mesh = Mesh()
    with HDF5File(MPI.comm_world, filename,'r') as f:
        f.read(mesh, "mesh", False)
        p_elem = eval(f.attributes("/pressure").to_dict()["signature"])
        v_elem = eval(f.attributes("/velocity").to_dict()["signature"])
        DG =     velocity_histo(sm, v, f"{dirname}/velocity_histo.png")
        FunctionSpace(mesh, p_elem)
        V = FunctionSpace(mesh, v_elem)
        DG0 = FunctionSpace(mesh, "DG", 0)
        p, v, sm = Function(DG), Function(V), Function(DG0)
        f.read(p, "pressure")
        f.read(v, "velocity")
        f.read(sm, "label")
    p = interpolate(p, FunctionSpace(mesh, "DG", v_elem.degree()))

    from IPython import embed; embed()
    velocity_histo(sm, v, f"{dirname}/velocity_histo.png")
    topology, cell_types, x = create_vtk_structures(V)
    grid = pv.UnstructuredGrid(topology, cell_types, x)
    grid["v"] = v.vector()[:].reshape(-1, 3)
    grid["p"] = p.vector()[:] - p.vector().min()
    if "cardiac" in dirname:
        grid["p"] *= (1 + alpha**2 / 8)
        config = cardiac_config
    else:
        config = prod_config

    gridsl = grid.slice(generate_triangles=True)
    streamlines = gridsl.streamlines(
        source_center=(0.0832,0.095, 0.076), 
        source_radius=0.08,
    max_step_length=0.5, vectors="v",
    surface_streamlines=True, interpolator_type="p", 
    compute_vorticity=False, max_time=10, n_points=40000,
    )

    p_bar_args=dict(title="p (Pa)", vertical=False,
                        height=0.08, width=0.6, position_x=0.2,
                        position_y=-0.0, title_font_size=52,
                        bold=False, font_family="times",
                        label_font_size=44, fmt="%.3g")
    v_bar_args = p_bar_args.copy()
    v_bar_args.update(dict(title="u (m/s)", vertical=True),
                            height=0.4, width=0.08, position_x=0.85,
                            position_y=0.08)

    pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
    pl.add_mesh(grid.clip(), scalars="p", clim=[0, np.round(grid["p"].max(),3)],
        cmap=config["pcmap"],
        scalar_bar_args=p_bar_args)
    vmax = np.round(np.linalg.norm(grid["v"], axis=1).max(), 3)
    pl.add_mesh(streamlines, line_width=2, cmap=config["vcmap"], clim=[0, vmax],
    scalar_bar_args=v_bar_args)
    pl.camera_position = 'yz'
    pl.camera.roll += 0
    pl.camera.azimuth += 30
    pl.camera.elevation += 10
    pl.camera.zoom(1.2)
    print(pl.camera.position)
    print(pl.camera.focal_point)
    pl.screenshot(f"{dirname}/csf_v.png",
    transparent_background=True)

if __name__ == "__main__":
    typer.run(plot_csf_flow)

