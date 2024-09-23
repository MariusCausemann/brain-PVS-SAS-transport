from fenics import *
import numpy as np
import typer
from pathlib import Path
from IPython import embed

def plot_R(R, filename):
    from vtk_adapter import create_vtk_structures
    import pyvista as pv
    from generate_synthseg_mesh import PARID

    V = R.function_space()
    topology, cell_types, x = create_vtk_structures(V)
    grid = pv.UnstructuredGrid(topology, cell_types, x)
    grid["R"] = R.vector()[:]

    bg = pv.read_meshio("mesh/standard/standard.xdmf")
    pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
    pl.add_mesh(bg.extract_cells(np.isin(bg["label"], [PARID])).clip(),
                color="black")
    pl.add_mesh(grid.clip(), log_scale=True, cmap="RdGy_r", 
                #clim=(0, R.vector().max()),
    scalar_bar_args=dict(title="R", vertical=False,
                         height=0.08, width=0.6, position_x=0.2,
                         position_y=-0.0, title_font_size=52,
                         bold=True, #font_family="arial",
                         label_font_size=44, fmt="%.1g"))
    pl.camera_position = 'yz'
    pl.camera.roll += 0
    pl.camera.azimuth += 30
    pl.camera.zoom(1.2)
    pl.screenshot(filename, transparent_background=True)

rho = 993 # kg/m^3
nu = 7e-7 # m^2/s
omega = 2*np.pi
h = 3e-3 / 2
alpha = np.sqrt(h**2 * omega / nu)

def get_dispersion_enhancement(csf_pressure_file:str, outfile:str):

    assume_unsteady = True
    mesh = Mesh()
    with XDMFFile(csf_pressure_file) as file:
        file.read(mesh)
        DG = FunctionSpace(mesh, "DG", 1)
        pressure_csf = Function(DG)
        file.read_checkpoint(pressure_csf, "pressure")

    #p_cont = project(pressure_csf, V)
    gradp = sqrt(inner(grad(pressure_csf), grad(pressure_csf)))
    P = gradp /(rho*omega*nu/h)
    if assume_unsteady:
        print(f"adjusting pressure by 1 + alpha**2 / 8 = {1 + alpha**2 / 8}")
        P *= (1 + alpha**2 / 8) # scale pressure with womersley
        R = P**2 / (alpha**3)
    else:
        R = P**2
    R_unsmoothed = project(R, DG)

    V = FunctionSpace(mesh, "CG", 1)
    u,v = TrialFunction(V), TestFunction(V)
    a = (Constant(1e-4)*inner(grad(u), grad(v)) + Constant(1)*u*v)*dx
    L = R*v*dx
    R = Function(V)
    solve(a==L, R)
    R = interpolate(R, DG)
    assert R.vector().min() > 0

    with XDMFFile(outfile) as file:
        file.write(mesh)
        file.write_checkpoint(R, "R", append=True)
        file.write_checkpoint(R_unsmoothed, "R_unsmoothed", append=True)

    plot_R(R, Path(outfile).with_suffix(".png"))
    #plot_R(R_unsmoothed, Path(outfile).with_suffix(".png"))


if __name__ == "__main__":
    typer.run(get_dispersion_enhancement)