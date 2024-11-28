from fenics import *
from xii import *
from plotting_utils import get_result_fenics
import typer
from typing import List
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config, minmax
from solver import pcws_constant
import numpy as np
import seaborn as sns
from IPython import embed
import pandas as pd
import yaml
from branch_marking import color_branches
from label_arteries import pointlabels
from test_map_on_global_coords_shift import map_kdtree
from peristalticflow import mesh_to_weighted_graph, nx
from generate_synthseg_mesh import CSFID, PARID, LVID, V34ID, CSFNOFLOWID
import joypy
from cmap import Colormap
from scipy.stats import binned_statistic

mol2mmol = 1e3

def compare_concentrations(modelname:str):
    config = read_config(f"configfiles/{modelname}.yml")
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, 3*dt*config["output_frequency"])

    data = read_config(f"results/{modelname}/mean_concentrations.yml")

    labels = ["CSF", "parenchyma", "PVS artery", "PVS vein"]
    tot_values = [[data[f"ctot_{l}_{t}"] for t in times] for l in labels]
    concentrations = [[data[f"cmean_{l}_{t}"] for t in times] for l in labels]
    labels = ["CSF", "PAR", "PVS artery", "PVS vein"]

    set_plotting_defaults()
    #sns.set_palette("BrBG", n_colors=4)
    sns.set_palette(["#0a9396","#94d2bd", "#e9d8a6", "#ee9b00"])

    fig, ax = plt.subplots(figsize=(5,4))
    tot = np.zeros_like(times)

    for q, l in zip(tot_values, labels):
        new_tot = tot + np.array(q)*mol2mmol
        fill = ax.fill_between(times/ (60*60), new_tot, tot, 
                               alpha=0.7, label=l)
        tot = new_tot
    
    plt.xlabel("time (h)")
    #ax2.set_ylabel('mean tracer concentration (mmol/l)')
    ax.set_ylabel("total tracer content (mmol)")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=4, 
              columnspacing=0.5, frameon=False)
    plt.tight_layout()
    filename = f"plots/{modelname}/{modelname}_total_conc.png"
    plt.savefig(filename, dpi=300)

    fig, ax = plt.subplots(figsize=(5,4))
    tot = np.zeros_like(times)

    for c, l in zip(concentrations, labels):
        new_tot = tot + np.array(q)*mol2mmol
        plt.plot(times / (60*60), c, "-", label=l, lw=5)
        tot = new_tot
    
    plt.xlabel("time (h)")
    #ax2.set_ylabel('mean tracer concentration (mmol/l)')
    ax.set_ylabel("mean tracer concentration (mmol/l)")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=4, 
              columnspacing=0.5, frameon=False)
    plt.tight_layout()
    filename = f"plots/{modelname}/{modelname}_mean_conc.png"
    plt.savefig(filename, dpi=300)




if __name__ == "__main__":
    typer.run(compare_concentrations)