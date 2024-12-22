import typer
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config
import numpy as np
import seaborn as sns
import os
from typing import List

mol2mmol = 1e3

def compare_concentrations(models:List[str]):

    config = read_config(f"configfiles/{models[0]}.yml")
    v = '_'.join(models)
    os.makedirs(f"plots/comparisons/{v}", exist_ok=True)

    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, 3*dt*config["output_frequency"])
    labels = ["CSF", "parenchyma", "PVS artery", "PVS vein"]
    results = dict()
    for m in models:
        print(f"reading model {m}")
        data = read_config(f"results/{m}/mean_concentrations.yml")
        tot_values = {l:[data[f"ctot_{l}_{t}"] for t in times] for l in labels}
        concentrations = {l:[data[f"cmean_{l}_{t}"] for t in times] for l in labels}
        results[m] = dict(conc=concentrations, tot=tot_values)

    set_plotting_defaults()
    #sns.set_palette("BrBG", n_colors=4)
    #sns.set_palette(["#0a9396","#94d2bd", "#e9d8a6", "#ee9b00"])
    colors = ["#0a9396","#81b8a5", "#ee6700", "#eea700"]
    coldict = {"modelA-strongVM":"#6610f2","modelA":"#0b774d", "modelA-PVS-disp":"#f7b801",
               "LargePVS":"#f7b801", "LargePVSA":"darkred"}
    namedict = {"modelA": "baseline", "modelA-strongVM":"high PVS flow",
                "LargePVS":"high PVS flow (dilated)", "LargePVSA":"baseline (dilated)"}
    for t in ["tot", "conc"]:
        fig, axes = plt.subplots(2,2,figsize=(9,4),sharex=True, sharey=True)
        for i, ax, l in zip(range(len(labels)), axes.flatten(), labels):
            for m in models:
                ax.plot(times/3600, results[m][t][l], label=namedict[m] if i==0 else None, color=coldict[m])
            ax.set_title(l)
        plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=len(models)*2,
                     columnspacing=0.3, frameon=False)
        fig.supxlabel("time (h)", fontsize="medium", y=0.06,x= 0.55 )
        fig.supylabel("mean concentration (mmol/l)", fontsize="medium", x=0.03, y=0.6)
        plt.tight_layout()
        filename = f"plots/comparisons/{v}/{v}_{t}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight',)
    
    # make grouped, stacked barplot of total mass
    timepoints = np.array([1,3, 6,12,24])
    timeidx = np.where(np.isin(times, timepoints*3600))[0]
    r = np.arange(len(timeidx))
    namedict = {"modelA": "baseline", "modelA-strongVM":"high PVS flow",
                "LargePVS":"high PVS flow (dilated)", "LargePVSA":"baseline (dilated)"}

    barWidth = 1/(len(models) + 1)
    off = (len(models) -1)*0.5*barWidth 
    fig, ax = plt.subplots(figsize=(9,4))
    for i,m in enumerate(models):
        bottom = np.zeros_like(r, dtype=np.float64)
        for l,c in zip(labels, colors):
            tot = np.array([results[m]["tot"][l][ti] for ti in timeidx])*mol2mmol
            plt.bar(r - off + i*barWidth, tot, color=c, edgecolor='white', 
                    width=barWidth, bottom=bottom, label=l if i==0 else None)
            bottom += tot
        plt.annotate(namedict[m], (r[0] - off + i*barWidth,bottom[0]), 
                     ha='center', va="bottom", rotation=90,
                      textcoords="offset points", xytext=(0, 5))
    plt.xticks(r, [f"{t}h" for t in timepoints])
    plt.ylabel("tracer content (mmol)")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=4, 
              columnspacing=0.5, frameon=False,  )
    filename = f"plots/comparisons/{v}/{v}_tot_bar.png"
    plt.savefig(filename, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    typer.run(compare_concentrations)