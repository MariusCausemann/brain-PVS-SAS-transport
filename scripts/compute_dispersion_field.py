from fenics import *
import numpy as np
import typer
from pathlib import Path
from vtk.util.numpy_support import numpy_to_vtk
from plot_subdomains import get_camera
from plotting_utils import set_plotting_defaults, read_config
from subdomain_ids import CSFID, LVID, V34ID
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import PercentFormatter
import seaborn as sns
import os 
import pandas as pd

def R_histo(sm, R, filename):
    set_plotting_defaults()
    sns.set_context('paper', font_scale=1.8)
    DG0 = sm.function_space()
    cellvol = assemble(TestFunction(DG0)*dx)
    domains = [[CSFID], [V34ID],]
    names = ["SAS", "V3 & V4"]
    indices = [np.isin(sm.vector()[:], dom) for dom in domains]      
    fig, axes = plt.subplots(ncols=len(domains), figsize=(7,3))
    nbins = 20
    colors = sns.color_palette("crest", n_colors=len(domains))
    
    # Initialize dictionary to collect histogram data across subplots
    source_data = {}

    for ax, name, idx, col in zip(axes, names, indices, colors):
        Rvals = R.vector()[idx]
        xlims = (0, np.percentile(Rvals, 80))
        bins_definition = np.logspace(start=np.log10(1e-1), stop=np.log10(200), num=nbins)
        
        # 1. Capture hist outputs (renamed variables to avoid clashing)
        counts, bin_edges, _ = ax.hist(
            Rvals, bins=bins_definition, weights=cellvol[idx] / cellvol[idx].sum(),
            color=col, alpha=0.7, label=name
        )
        
        # 2. Calculate bin centers for the X-axis mapping
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        
        ax.set_xscale("log")
        mean = np.average(Rvals, weights=cellvol[idx])
        ax.axvline(mean, 
                   color="black", ymax=0.6, label=f"{mean:.1f} \n (mean)")
        
        ax.set_xlabel(f"R")
        ax.set_ylabel("frequency")
        ax.set_xticks([0.1, 1, 10, 100])
        ax.set_ylim((0, max(counts)*1.5))
        ax.legend(frameon=False, loc="lower left", alignment="right",
                  borderaxespad=0, handlelength=1.3, bbox_to_anchor=[0.2, 0.7])
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        
        # 3. Store the data for this specific subplot panel
        safe_name = name.replace(" & ", "_")
        source_data[f"{safe_name}_bin_center"] = bin_centers
        source_data[f"{safe_name}_relative_frequency"] = counts
        
        # Pad the single scalar mean value with NaNs to match column lengths
        mean_column = np.full(len(counts), np.nan)
        mean_column[0] = mean
        source_data[f"{safe_name}_mean_R"] = mean_column

    # Save original plot graphic
    plt.savefig(filename, bbox_inches='tight')
    
    # 4. Generate and save the accompanying Source Data CSV
    base_path, _ = os.path.splitext(filename)
    csv_filename = f"{base_path}_source_data.csv"
    
    df = pd.DataFrame(source_data)
    df.to_csv(csv_filename, index=False)
    print(f"--> Exported histogram source data to: {csv_filename}")


def plot_R(R, filename, background_mesh_file):
    from vtk_adapter import create_vtk_structures
    import pyvista as pv
    from subdomain_ids import PARID

    V = R.function_space()
    topology, cell_types, x = create_vtk_structures(V)
    grid = pv.UnstructuredGrid(topology, cell_types, x)
    grid["R"] = R.vector()[:]

    Rmin, Rmax = grid["R"].min(), grid["R"].max()
    bg = pv.read_meshio(background_mesh_file)
    if "respiratory" in str(filename):
        cmap = sns.blend_palette(["black", "white", "#092f8d"],as_cmap=True)
    elif "cardiac" in str(filename):
        cmap = "RdGy_r"

    for c in ["white", "black"]:

        pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
        pl.add_mesh(bg.extract_cells(np.isin(bg["label"], [PARID])).clip(),
                color="black")

        pl.add_mesh(grid.clip(), log_scale=True,cmap=cmap,
                    show_scalar_bar=False, clim=(Rmin, Rmax))
        bar = pl.add_scalar_bar(title="R", vertical=False,
                        height=0.08, width=0.6, position_x=0.2,
                        position_y=-0.0, title_font_size=64,
                        bold=False, font_family="times", color=c,
                        label_font_size=60, fmt="%.1f")
        ticks = [Rmin, 1e-1,1, 10, int(np.floor(Rmax / 10) * 10)]
        bar.SetCustomLabels(numpy_to_vtk(ticks))
        bar.SetUseCustomLabels(True)
        pl.camera = get_camera(bg)
        path = Path(filename)
        if c=="white":
            path = path.parent / (path.stem + "_dark" + path.suffix)
        pl.screenshot(path, transparent_background=True)


rho = 993 # kg/m^3
nu = 7e-7 # m^2/s
h = 3e-3 / 2

omega_cardiac = 2*np.pi * 1.0
omega_respiratory = 2*np.pi * 0.25
alpha_cardiac = np.sqrt(h**2 * omega_cardiac / nu)
alpha_respiratory = np.sqrt(h**2 * omega_respiratory / nu)

def test_sharp2019_SSAS():
    delta_p = 45.3 # Pa
    l = 0.1 # m
    h = 3e-3 / 2
    omega = omega_cardiac
    gradp = delta_p / l
    P = gradp /(rho*omega*nu/h)
    assert np.isclose(P, 155.7, rtol=5e-2) # value from paper
    alpha = np.sqrt(h**2 * omega_cardiac / nu)
    R = P**2 / (alpha**3)

def PVS_R():
    delta_p = 2.5 # Pa
    l = 0.1 # m
    omega = 2*np.pi
    gradp = delta_p / l
    h = 1e-3 / 2
    P = gradp /(rho*omega*nu/h)
    alpha = np.sqrt(h**2 * omega / nu)
    R = P**2 / (alpha**3)


def get_dispersion_enhancement(csf_pressure_file:str, outfile:str):
    if "cardiac" in csf_pressure_file:
        omega = omega_cardiac
        alpha = alpha_cardiac
    elif "respiratory" in csf_pressure_file:
        omega = omega_respiratory
        alpha = alpha_respiratory
    assume_unsteady = True
    mesh = Mesh()
    with HDF5File(MPI.comm_world, csf_pressure_file,'r') as f:
        f.read(mesh, "mesh", False)
        p_elem = eval(f.attributes("/pressure").to_dict()["signature"])
        DG = FunctionSpace(mesh, p_elem)
        pressure_csf = Function(DG)
        f.read(pressure_csf, "pressure")
        DG0 = FunctionSpace(mesh, "DG", 0)
        sm = Function(DG0)
        f.read(sm, "label")

    gradp = sqrt(inner(grad(pressure_csf), grad(pressure_csf)))
    P = gradp /(rho*omega*nu/h)
    print(f"alpha**2: {alpha**2}")

    if assume_unsteady:
        print(f"adjusting pressure by 1 + alpha**2 / 8 = {1 + alpha**2 / 8}")
        P *= (1 + alpha**2 / 8) # scale pressure with womersley
        R = P**2 / (alpha**3)
    else:
        R = P**2
    R_unsmoothed = project(R, DG, solver_type="cg", preconditioner_type="petsc_amg")

    V = FunctionSpace(mesh, "CG", 1)
    u,v = TrialFunction(V), TestFunction(V)
    a = (Constant(1e-4)*inner(grad(u), grad(v)) + Constant(1)*u*v)*dx
    L = R*v*dx
    R = Function(V)
    solve(a==L, R, solver_parameters={"linear_solver":"cg", "preconditioner":"petsc_amg"})
    R = interpolate(R, FunctionSpace(mesh, "DG", 1))
    if R.vector().min() < 0:
        print(f"warning: R_min = {R.vector().min()}")
        R.vector()[:] -= R.vector().min()
    print(f"R max: {R.vector().max()}")

    R_histo(sm, interpolate(R, DG0), Path(outfile).parent / "R_histo.svg")

    with XDMFFile(outfile) as file:
        file.write(mesh)
        file.write_checkpoint(R, "R", append=True)
        file.write_checkpoint(R_unsmoothed, "R_unsmoothed", append=True)

    config = read_config(f"configfiles/{Path(csf_pressure_file).parent.name}.yml")
    meshname = config["mesh"].removesuffix("_outer.xdmf")
    background_mesh_file = f"{meshname}.xdmf"
    plot_R(R, Path(outfile).with_suffix(".png"), background_mesh_file)


if __name__ == "__main__":
    typer.run(get_dispersion_enhancement)