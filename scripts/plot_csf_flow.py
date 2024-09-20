
import pyvista as pv
import typer
from vtk_adapter import create_vtk_structures
from solver  import get_mesh, read_vtk_network, as_P0_function
from pathlib import Path
import numpy as np
from fenics import *
import matplotlib 
import k3d.colormaps.paraview_color_maps as pcm
from compute_dispersion_field import alpha

def from_k3d(colorlist):
    cm = np.array(colorlist).reshape(-1, 4)[:, 1:]
    return matplotlib.colors.LinearSegmentedColormap.from_list("",cm)


cardiac_config = dict(pcmap="balance", vcmap="turbo")
prod_config = dict(pcmap=from_k3d([1]*12 + pcm.Blue___Green___Orange), 
vcmap="inferno")

def plot_csf_flow(dirname: str):
    filename_v = f"{dirname}/csf_vis_v.xdmf"
    filename_p = f"{dirname}/csf_vis_p.xdmf"

    mesh = Mesh()
    with XDMFFile(filename_p) as f:
        f.read(mesh)
        DG = FunctionSpace(mesh, "DG", 1)
        p = Function(DG)
        f.read_checkpoint(p, "pressure")
        p = interpolate(p, FunctionSpace(mesh, "DG", 2))

    with XDMFFile(filename_v) as f:
        V = VectorFunctionSpace(mesh, "DG", 2)
        v = Function(V)
        f.read_checkpoint(v, "velocity")

    topology, cell_types, x = create_vtk_structures(V)
    grid = pv.UnstructuredGrid(topology, cell_types, x)
    grid["v"] = v.vector()[:].reshape(-1, 3)
    grid["p"] = p.vector()[:]
    if "cardiac" in dirname:
        grid["p"] *= alpha**2
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
                        bold=True, #font_family="arial",
                        label_font_size=44, fmt="%.3g")
    v_bar_args = p_bar_args.copy()
    v_bar_args.update(dict(title="u (m/s)", vertical=True),
                            height=0.4, width=0.08, position_x=0.85,
                            position_y=0.08)

    pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
    pl.add_mesh(grid.clip(), scalars="p", clim=[0, np.round(grid["p"].max(),2)],
        cmap=config["pcmap"],
        scalar_bar_args=p_bar_args)
    vmax = np.round(np.linalg.norm(grid["v"], axis=1).max(), 3)
    pl.add_mesh(streamlines, line_width=2, cmap=config["vcmap"], clim=[0, vmax],
    scalar_bar_args=v_bar_args)
    pl.camera_position = 'yz'
    pl.camera.roll += 0
    pl.camera.azimuth += 30
    pl.camera.zoom(1.2)
    pl.screenshot(f"{dirname}/csf_v.png",
    transparent_background=True)

if __name__ == "__main__":
    typer.run(plot_csf_flow)

