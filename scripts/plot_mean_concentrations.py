from fenics import *
import typer
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config, minmax
import numpy as np
import seaborn as sns
import os

mol2mmol = 1e3

def compare_concentrations(modelname:str):
    config = read_config(f"configfiles/{modelname}.yml")
    os.makedirs(f"plots/{modelname}", exist_ok=True)
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, dt*config["output_frequency"])

    data = read_config(f"results/{modelname}/mean_concentrations.yml")

    labels = ["CSF", "parenchyma", "PVS artery", "PVS vein",]
    if f"ctot_Tips_0" in data.keys(): labels.append("Tips")

    tot_values = [[data[f"ctot_{l}_{t}"] for t in times] for l in labels]
    concentrations = [[data[f"cmean_{l}_{t}"] for t in times] for l in labels]
    labels = ["CSF", "PAR", "PVS art.", "PVS vein", "PVS cont."]

    set_plotting_defaults()
    #sns.set_palette("BrBG", n_colors=4)
    #sns.set_palette(["#0a9396","#94d2bd", "#e9d8a6", "#ee9b00"])
    sns.set_palette(["#0a9396","#81b8a5", "#ee6700", "#f3e30d","#eea700",])

    fig, ax = plt.subplots(figsize=(4.5,3.5))
    tot = np.zeros_like(times)
    prev = -100
    for q, l in zip(tot_values, labels):
        new_tot = tot + np.array(q)*mol2mmol
        fill = ax.fill_between(times/ (60*60), new_tot, tot, 
                               alpha=1, label=l)
        ypos =  max(tot[-1], prev+0.04)
        final_tot = sum(t[-1] for t in tot_values)*mol2mmol
        #plt.scatter(times[-1]/(60*60), ypos, marker="o", color=fill.get_facecolor())
        plt.annotate(f"{100*mol2mmol*q[-1]/final_tot:.0f}%",(times[-1]/(60*60), ypos),
                     color=fill.get_facecolor(),
                     xytext=(2, 0), textcoords='offset points',
                     xycoords="data",)
        print(f"m_final {l}: {q[-1]}")
        prev = ypos
        tot = new_tot

    print(tot[-5:].min())    
    print(tot[-5:].max())    
    plt.xlabel("time (h)")
    #ax2.set_ylabel('mean tracer concentration (mmol/l)')
    ax.set_ylabel("total tracer content (mmol)")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3, 
              columnspacing=0.5, frameon=False)
    plt.tight_layout()
    filename = f"plots/{modelname}/{modelname}_total_conc.png"
    plt.savefig(filename, dpi=300, transparent=True)

    fig, ax = plt.subplots(figsize=(4.5,3.5))
    tot = np.zeros_like(times)

    for c, l in zip(concentrations, labels):
        new_tot = tot + np.array(q)*mol2mmol
        plt.plot(times / (60*60), c, "-", label=l, lw=5)
        tot = new_tot
        tmax = times[np.argmax(c)]
        print(f"c_max: {l}: {max(c)} at {tmax / 60} min")
    
    plt.xlabel("time (h)")
    #ax2.set_ylabel('mean tracer concentration (mmol/l)')
    ax.set_ylabel("mean tracer concentration (mmol/l)")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3, 
              columnspacing=0.5, frameon=False)
    plt.tight_layout()
    filename = f"plots/{modelname}/{modelname}_mean_conc.png"
    plt.savefig(filename, dpi=300, transparent=True)




if __name__ == "__main__":
    typer.run(compare_concentrations)