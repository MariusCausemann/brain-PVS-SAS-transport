
from xml import dom
import pyvista as pv
import typer
from vtk_adapter import create_vtk_structures
import numpy as np
from fenics import *
import matplotlib 
import k3d.colormaps.paraview_color_maps as pcm
from compute_dispersion_field import alpha_cardiac, alpha_respiratory
from plotting_utils import set_plotting_defaults, read_config
from subdomain_ids import CSFID, LVID, V34ID
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import seaborn as sns
from plot_subdomains import get_camera
from vtk.util.numpy_support import numpy_to_vtk
import seaborn as sns
from pathlib import Path
from mark_and_refine_mesh import SPINAL_OUTLET_ID

m2mum = 1e6

def velocity_histo(sm, u, filename, vscale, vunit):
    set_plotting_defaults()
    sns.set_context('paper', font_scale=1.8)

    DG0 = sm.function_space()
    umag = project(sqrt(inner(u,u)), DG0,
                   solver_type="cg", preconditioner_type="hypre_amg")
    cellvol = assemble(TestFunction(DG0)*dx)
    domains = [[CSFID], [V34ID],]
    names = ["SAS", "V3 & V4"]
    indices = [np.isin(sm.vector()[:], dom) for dom in domains]  
    #from IPython import embed; embed()
    fig, axes = plt.subplots(ncols=len(domains), figsize=(7,3))
    nbins = 20
    colors = sns.color_palette("crest", n_colors=len(domains))
    for ax, name, idx, col in zip(axes,names , indices, colors):
        uvals = umag.vector()[idx] * vscale
        xlims = (0, np.percentile(uvals, 92))
        ax.hist(uvals, bins=nbins, weights=cellvol[idx] / cellvol[idx].sum(),
                 range=xlims,
                  color=col, alpha=0.7, label=name)
        mean = np.average(uvals, weights=cellvol[idx])
        ax.axvline(mean, 
                   color="black",ymax=0.6, label=f"{mean:.1f} {vunit} \n (mean)")
        #ax.axvline(np.max(uvals), ls="dotted",
        #           color="black",ymax=0.6, label=f"{max*1e-3:.1f} mm/s")
        ax.set_xlabel(f"velocity ({vunit})")
        ax.set_ylabel("frequency")
        #plt.tight_layout()
        ax.set_xlim(xlims)
        ax.legend(frameon=False, loc="lower left", alignment="right",
                  borderaxespad=0, handlelength=1.3, bbox_to_anchor=[0.1, 0.7])
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    plt.savefig(filename, bbox_inches='tight')

def from_k3d(colorlist):
    cm = np.array(colorlist).reshape(-1, 4)[:, 1:]
    return matplotlib.colors.LinearSegmentedColormap.from_list("", cm)

cardiac_config = dict(pcmap="balance", 
                    vcmap="inferno", vscale=1e2, pscale=1, 
                    vunit="cm/s", punit="Pa",
                    vticks= [1e-2, 0.1, 1, 5],
                    vmin=1e-2, vbar_fmt="%.1f")
respiration_config = dict(pcmap="balance", 
                    vcmap="inferno", vscale=1e2, pscale=1, 
                    vunit="cm/s", punit="Pa",
                    vticks= [1e-2, 0.1, 1, 5],
                    vmin=1e-2, vbar_fmt="%.1f")
prod_config = dict(pcmap=from_k3d([1]*4 + pcm.Blue___Green___Orange), 
                    vcmap="inferno", vscale=1e3, pscale=1e3,  
                    vunit="mm/s", punit="mPa", vticks= [1e-4, 1e-3, 1e-2, 1e-1],
                    vmin=1e-4,  vbar_fmt="%.3f")

def plot_csf_flow(dirname: str):
    filename = f"{dirname}/flow.hdf"

    mesh = Mesh()
    with HDF5File(MPI.comm_world, filename,'r') as f:
        f.read(mesh, "mesh", False)
        p_elem = eval(f.attributes("/pressure").to_dict()["signature"])
        v_elem = eval(f.attributes("/velocity").to_dict()["signature"])
        DG = FunctionSpace(mesh, p_elem)
        V = FunctionSpace(mesh, v_elem)
        DG0 = FunctionSpace(mesh, "DG", 0)
        p, v, sm = Function(DG), Function(V), Function(DG0)
        f.read(p, "pressure")
        f.read(v, "velocity")
        f.read(sm, "label")
    p = interpolate(p, FunctionSpace(mesh, "DG", v_elem.degree()))

    flow_model = Path(dirname).name
    flow_config = read_config(f"configfiles/{flow_model}.yml")
    bm = MeshFunction("size_t", mesh, 2, 0)
    with XDMFFile(flow_config["mesh"].replace(".xdmf","_facets.xdmf")) as f:
        f.read(bm, "f")

    ds = Measure("ds", mesh, subdomain_data=bm)
    n = FacetNormal(mesh)
    fm_area = assemble(Constant(1)*ds(SPINAL_OUTLET_ID))
    fm_outflow = assemble(inner(v, n)*ds(SPINAL_OUTLET_ID))
    fm_mean_velocity = fm_outflow/ fm_area
    print(fm_mean_velocity)
    topology, cell_types, x = create_vtk_structures(V)
    grid = pv.UnstructuredGrid(topology, cell_types, x)
    grid["v"] = v.vector()[:].reshape(-1, 3)
    grid["p"] = p.vector()[:] - p.vector().min()
    if "cardiac" in dirname:
        print("adjusting pressure with (1 + alpha_cardiac**2 / 8) ")
        print(f"alpha_cardiac = {alpha_cardiac} ")
        grid["p"] *= (1 + alpha_cardiac**2 / 8)
        config = cardiac_config
        velocity_histo(sm, v, f"{dirname}/cardiac_velocity_histo.png", vscale=1e3,
                       vunit="mm/s")
    elif "respiratory" in dirname:
        print("adjusting pressure with (1 + alpha_respiratory**2 / 8) ")
        print(f"alpha_respiratory = {alpha_respiratory} ")
        grid["p"] *= (1 + alpha_respiratory**2 / 8)
        config = respiration_config
        velocity_histo(sm, v, f"{dirname}/resp_velocity_histo.png", vscale=1e3,
                       vunit="mm/s")
    else:
        config = prod_config
        velocity_histo(sm, v, f"{dirname}/prod_velocity_histo.png", vscale=1e6,
                       vunit="μm/s")

    grid["v"] *= config["vscale"]
    grid["p"] *= config["pscale"]

    gridsl = grid.slice(generate_triangles=True)
    streamlines = gridsl.streamlines(
        source_center=(0.0832,0.095, 0.076), 
        source_radius=0.08,
    max_step_length=0.5, vectors="v",
    surface_streamlines=True, interpolator_type="p", 
    compute_vorticity=False, n_points=40000, 
    terminal_speed=3e-4, max_steps=200,
    )

    p_bar_args=dict(title=f" p ({config['punit']})", vertical=False,
                        height=0.08, width=0.6, position_x=0.2,
                        position_y=-0.0, title_font_size=64,
                        bold=False, font_family="times",
                        label_font_size=60, fmt="%.1f")
    v_bar_args = p_bar_args.copy()
    v_bar_args.update(dict(title=f"u ({config['vunit']})", vertical=True),
                            height=0.4, width=0.08, position_x=0.87,
                            position_y=0.02,
                            fmt=config["vbar_fmt"])
    for c in ["white", "black"]:

        pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
        pl.add_mesh(grid.clip(), scalars="p", clim=[0, np.round(grid["p"].max(),1)],
            cmap=config["pcmap"],
            scalar_bar_args=dict(**p_bar_args, color=c))
        vmax = np.round(np.linalg.norm(grid["v"], axis=1).max(), 3)
        pl.add_mesh(streamlines, line_width=2, cmap=config["vcmap"],
                    clim=[config["vmin"], vmax],show_scalar_bar=False,
                    log_scale=True)
        
        bar = pl.add_scalar_bar(**v_bar_args, color=c)
        ticks = config["vticks"] + [np.round(vmax, 3)]
        bar.SetCustomLabels(numpy_to_vtk(ticks))
        bar.SetUseCustomLabels(True)
        pl.camera_position = 'yz'
        pl.camera = get_camera(grid)
        pl.screenshot(f"{dirname}/csf_v{'_dark' if c=='white' else ''}.png",
        transparent_background=True)

    #from IPython import embed; embed()

if __name__ == "__main__":
    typer.run(plot_csf_flow)

