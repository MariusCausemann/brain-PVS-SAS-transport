import typer
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config, minmax
import numpy as np
import pyvista as pv

mol2mmol = 1e3

def plot_tip_conc(modelname:str):
    config = read_config(f"configfiles/{modelname}.yml")
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, dt*config["output_frequency"])
    tip_conc = pv.read(f"results/{modelname}/{modelname}_tips.vtk")

    c_vals = np.array([tip_conc[f"c_tips_{t}"] for t in times])

    set_plotting_defaults()
    plt.figure(figsize=(4.3,3.2))
    ivd_lines = plt.plot(times / (60*60), c_vals, alpha=0.2, color="black", 
             linewidth=0.6)
    mean_line = plt.plot(times / (60*60), c_vals.mean(axis=1), color="orange")
    plt.xlabel("time (h)")
    plt.ylabel("concentration (mmol/l)")
    plt.legend([ivd_lines[0], mean_line[0]], 
           ["0D compartments", "Average concentration"], frameon=False)
    plt.tight_layout()
    plt.savefig(f"plots/{modelname}/{modelname}_tip_conc.png", dpi=300,
                transparent=True)

def plot_leaf_conc(modelname:str):
    plot_tip_conc(modelname)
    config = read_config(f"configfiles/{modelname}.yml")
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, dt*config["output_frequency"])
    import dolfin as df
    from plotting_utils import get_result_fenics
    art_conc, art_radii, art_roots = get_result_fenics(modelname, "artery", times, getroots=True)

    tip_ids = np.where(art_roots.array()==1)[0]
    tip_dofs = df.vertex_to_dof_map(art_conc[0].function_space())[tip_ids]
    c_vals = np.array([c.vector()[tip_dofs] for c in art_conc])

    set_plotting_defaults()
    plt.figure(figsize=(4.3,3.2))
    ivd_lines = plt.plot(times / (60*60), c_vals, alpha=0.2, color="black", 
             linewidth=0.6)
    mean_line = plt.plot(times / (60*60), c_vals.mean(axis=1), color="red")
    plt.xlabel("time (h)")
    plt.ylabel("concentration (mmol/l)")
    plt.legend([ivd_lines[0], mean_line[0]], 
           ["leaf concentration", "Average concentration"], frameon=False)
    plt.tight_layout()
    plt.savefig(f"plots/{modelname}/{modelname}_leaf_conc.png",
                 dpi=300, transparent=True)

if __name__ == "__main__":
    typer.run(plot_leaf_conc)