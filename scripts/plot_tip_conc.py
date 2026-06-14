import typer
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config, minmax
import numpy as np
import pyvista as pv
import os
import pandas as pd

mol2mmol = 1e3

def plot_tip_conc(modelname:str):
    config = read_config(f"configfiles/{modelname}.yml")
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, dt*config["output_frequency"])
    time_hours = times / (60 * 60)
    
    tip_conc = pv.read(f"results/{modelname}/{modelname}_tips.vtk")
    c_vals = np.array([tip_conc[f"c_tips_{t}"] for t in times])

    set_plotting_defaults()
    plt.figure(figsize=(4.3,3.2))
    ivd_lines = plt.plot(time_hours, c_vals, alpha=0.2, color="black", 
             linewidth=0.6)
    mean_curve = c_vals.mean(axis=1)
    mean_line = plt.plot(time_hours, mean_curve, color="orange")
    
    plt.xlabel("time (h)")
    plt.ylabel("concentration (mmol/l)")
    plt.legend([ivd_lines[0], mean_line[0]], 
           ["0D compartments", "Average concentration"], frameon=False)
    plt.tight_layout()
    
    # Save Image
    img_filename = f"plots/{modelname}/{modelname}_tip_conc.svg"
    plt.savefig(img_filename, dpi=300, transparent=True)
    plt.close()

    # ==========================================================================
    # GENERATE SOURCE DATA (0D Compartments)
    # ==========================================================================
    source_data = {"Time (h)": time_hours}
    
    # Dynamically unpack every individual compartment trace into its own column
    num_compartments = c_vals.shape[1]
    for idx in range(num_compartments):
        source_data[f"Compartment_{idx+1}_concentration_mmol_l"] = c_vals[:, idx]
        
    # Add the final orange average trace
    source_data["Average_concentration_mmol_l"] = mean_curve
    
    base_path, _ = os.path.splitext(img_filename)
    pd.DataFrame(source_data).to_csv(f"{base_path}_source_data.csv", index=False)
    print(f"--> Exported 0D tip trajectories to: {base_path}_source_data.csv")


def plot_leaf_conc(modelname:str):
    plot_tip_conc(modelname)
    config = read_config(f"configfiles/{modelname}.yml")
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, dt*config["output_frequency"])
    time_hours = times / (60 * 60)
    
    import dolfin as df
    from plotting_utils import get_result_fenics
    art_conc, art_radii, art_roots = get_result_fenics(modelname, "artery", times, getroots=True)

    tip_ids = np.where(art_roots.array()==1)[0]
    tip_dofs = df.vertex_to_dof_map(art_conc[0].function_space())[tip_ids]
    c_vals = np.array([c.vector()[tip_dofs] for c in art_conc])

    set_plotting_defaults()
    plt.figure(figsize=(4.3,3.2))
    ivd_lines = plt.plot(time_hours, c_vals, alpha=0.2, color="black", 
             linewidth=0.6)
    mean_curve = c_vals.mean(axis=1)
    mean_line = plt.plot(time_hours, mean_curve, color="red")
    
    plt.xlabel("time (h)")
    plt.ylabel("concentration (mmol/l)")
    plt.legend([ivd_lines[0], mean_line[0]], 
           ["leaf concentration", "Average concentration"], frameon=False)
    plt.tight_layout()
    
    # Save Image (Fixed extension typo from .scg to .svg)
    img_filename = f"plots/{modelname}/{modelname}_leaf_conc.svg"
    plt.savefig(img_filename, dpi=300, transparent=True)
    plt.close()

    # ==========================================================================
    # GENERATE SOURCE DATA (Leaf Concentration)
    # ==========================================================================
    source_data = {"Time (h)": time_hours}
    
    # Dynamically unpack every individual leaf trace into its own column
    num_leaves = c_vals.shape[1]
    for idx in range(num_leaves):
        source_data[f"Leaf_{idx+1}_concentration_mmol_l"] = c_vals[:, idx]
        
    # Add the final red average trace
    source_data["Average_concentration_mmol_l"] = mean_curve
    
    base_path, _ = os.path.splitext(img_filename)
    pd.DataFrame(source_data).to_csv(f"{base_path}_source_data.csv", index=False)
    print(f"--> Exported leaf node trajectories to: {base_path}_source_data.csv")


if __name__ == "__main__":
    typer.run(plot_tip_conc)
    typer.run(plot_leaf_conc)
