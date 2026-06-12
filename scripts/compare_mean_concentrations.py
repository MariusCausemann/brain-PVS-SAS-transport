import typer
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config
import numpy as np
import seaborn as sns
import os
from typing import List
import pandas as pd

mol2mmol = 1e3

def compare_concentrations(models:List[str]):

    config = read_config(f"configfiles/{models[0]}.yml")
    v = '_'.join(models)
    os.makedirs(f"plots/comparisons/{v}", exist_ok=True)

    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, dt*config["output_frequency"])
    time_hours = times / 3600  # Baseline tracking timeline
    labels = ["CSF", "parenchyma", "PVS artery", "PVS vein"]
    results = dict()
    for m in models:
        print(f"reading model {m}")
        data = read_config(f"results/{m}/mean_concentrations.yml")
        tot_values = {l:[data[f"ctot_{l}_{t}"] for t in times] for l in labels}
        concentrations = {l:[data[f"cmean_{l}_{t}"] for t in times] for l in labels}
        results[m] = dict(conc=concentrations, tot=tot_values)

    set_plotting_defaults()
    colors = ["#0a9396","#81b8a5", "#ee6700", "#eea700"]
    coldict = {"modelA-strongVM":"#6610f2","modelA":"#0b774d", "modelA-PVS-disp":"#f7b801",
               "LargePVS":"#f7b801", "LargePVSA":"darkred"}
    namedict = {"modelA": "baseline", "modelA-strongVM":"high PVS flow",
                "LargePVS":"high PVS flow (dilated)", "LargePVSA":"baseline (dilated)"}
    
    # ==========================================================================
    # FIGURES 1 & 2: 2x2 Grids for Total Mass and Concentrations
    # ==========================================================================
    for t in ["tot", "conc"]:
        fig, axes = plt.subplots(2,2,figsize=(9,4),sharex=True, sharey=True)
        
        # Initialize dictionary to collect 2x2 subpanel vectors
        source_data_grid = {"Time (h)": time_hours}

        for i, ax, l in zip(range(len(labels)), axes.flatten(), labels):
            for m in models:
                ax.plot(time_hours, results[m][t][l], label=namedict[m] if i==0 else None, color=coldict[m])
                
                # Format clean tracking header strings (e.g., parenchyma_modelA_conc)
                safe_label = l.replace(" ", "_")
                metric_unit = "mmol" if t == "tot" else "mmol_l"
                source_data_grid[f"{safe_label}_{m}_{t}_{metric_unit}"] = results[m][t][l]

            ax.set_title(l)
            
        plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=len(models)*2,
                     columnspacing=0.3, frameon=False)
        fig.supxlabel("time (h)", fontsize="medium", y=0.06,x= 0.55 )
        fig.supylabel("mean concentration (mmol/l)" if t == "conc" else "total content", fontsize="medium", x=0.03, y=0.6)
        plt.tight_layout()
        
        # Save Image
        filename = f"plots/comparisons/{v}/{v}_{t}.svg"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save Grid Source Data CSV
        csv_filename = filename.replace(".svg", "_source_data.csv")
        pd.DataFrame(source_data_grid).to_csv(csv_filename, index=False)
        print(f"--> Exported time-series grid source data to: {csv_filename}")
    
    # ==========================================================================
    # FIGURE 3: Grouped, Stacked Bar Plot of Total Mass
    # ==========================================================================
    timepoints = np.array([1,3, 6,12,24])
    timeidx = np.where(np.isin(times, timepoints*3600))[0]
    r = np.arange(len(timeidx))

    barWidth = 1/(len(models) + 1)
    off = (len(models) -1)*0.5*barWidth 
    fig, ax = plt.subplots(figsize=(9,4))
    
    # Initialize dictionary to track stacked bar values at defined intervals
    source_data_bar = {"Timepoint": [f"{t}h" for t in timepoints]}

    for i,m in enumerate(models):
        bottom = np.zeros_like(r, dtype=np.float64)
        for l,c in zip(labels, colors):
            tot = np.array([results[m]["tot"][l][ti] for ti in timeidx])*mol2mmol
            plt.bar(r - off + i*barWidth, tot, color=c, edgecolor='white', 
                    width=barWidth, bottom=bottom, label=l if i==0 else None)
            
            # Map out each structural compartment layer thickness (magnitude) per model
            safe_label = l.replace(" ", "_")
            source_data_bar[f"{m}_{safe_label}_content_mmol"] = tot
            
            bottom += tot
            
        plt.annotate(namedict[m], (r[0] - off + i*barWidth,bottom[0]), 
                     ha='center', va="bottom", rotation=90,
                      textcoords="offset points", xytext=(0, 5))
                      
    plt.xticks(r, [f"{t}h" for t in timepoints])
    plt.ylabel("tracer content (mmol)")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=4, 
              columnspacing=0.5, frameon=False)
              
    # Save Image
    filename = f"plots/comparisons/{v}/{v}_tot_bar.svg"
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Save Bar Chart Source Data CSV
    csv_filename = filename.replace(".svg", "_source_data.csv")
    pd.DataFrame(source_data_bar).to_csv(csv_filename, index=False)
    print(f"--> Exported stacked bar source data to: {csv_filename}")


if __name__ == "__main__":
    typer.run(compare_concentrations)