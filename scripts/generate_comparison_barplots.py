import yaml
from plotting_utils import read_config
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
from plotting_utils import time_str, set_plotting_defaults
from compute_dispersion_field import alpha_cardiac, alpha_respiratory
import typer
from typing import List


mesh_refinement_models = ["standard", "vasomotion"]
time_refinement_models = ["dt=240s", "dt=120s", "dt=60s"]
variants = ["mesh_refinement", "time_refinement"]

def compare_models(models: List[str]):
    v = "_".join(models)
    os.makedirs(f"plots/comparisons/{v}/", exist_ok=True)

    results = dict()
    for m in models:
        md = dict()
        tm = m
        tmconfig = read_config(f"configfiles/{tm}.yml")
        actual_mesh = Path(tmconfig["mesh"])
        print(actual_mesh)
        # get mesh stats
        md.update(read_config(actual_mesh.with_suffix(".yml")))
        # get CSF flow stats
        csfm = tmconfig["csf_velocity_file"].split("/")[-2]
        md.update(read_config(f"results/csf_flow/{csfm}/metrics.yml"))
        # get cardiac PVS flow stats
        dispersion_file = tmconfig["csf_dispersion_file"]
        if isinstance(dispersion_file, list): dispersion_file=dispersion_file[0]
        csfm_cardiac = dispersion_file.split("/")[-2]
        cardiac_csf = read_config(f"results/csf_flow/{csfm_cardiac}/metrics.yml")
        cardiac_csf["pmax"] *= (1 + alpha_cardiac ** 2 / 8)
        md.update({f"cardiac_{k}":v for k,v in cardiac_csf.items()})
        # get mean concentrations
        md.update(read_config(f"results/{tm}/mean_concentrations.yml"))
        trmetrics = read_config(f"results/{tm}/{tm}_metrics.yml")
        md.update({"Rmax":trmetrics["R_max"],"Rmean":trmetrics["R_mean"],
                "Rmin":trmetrics["R_min"] })
        for dom in ["csf", "art", "ven"]:
            md.update({f"{dom}_max": max(trmetrics[f"{dom}_max"]),
                    f"{dom}_min": min(trmetrics[f"{dom}_min"]),
                    f"rel_undershot_{dom}": abs(min(np.array(trmetrics[f"{dom}_min"]) /
                                                np.array(trmetrics[f"{dom}_max"]) )) },
                    )
        results[m] = md

    df = pd.DataFrame(results)
    set_plotting_defaults()

    barplot_groups = [
    {"vars":["hmin", "hmax"],"varnames":[r"$\rm h_{min}$", r"$\rm h_{max}$"], 
    "ylabel":"mm", "title":"cell size", "scale":1e3},
    {"vars":["ncells", "npoints"],"varnames":["#cells", "#points"],
    "ylabel":" # millon", "title":"number of cells", "scale":1e-6},
    {"vars":["Rmean"],"varnames":[r"$\rm R_{mean}$"],
    "ylabel":"R", "title":"mean dispersion"},
    {"vars":["Rmax"],"varnames":[r"$\rm R_{max}$"],
    "ylabel":"R", "title":"max dispersion"},
    {"vars":["pmax"],"varnames":[r"$\rm p_{max}$"],
    "ylabel":"Pa", "title":"CSF pressure"},
    {"vars":["cardiac_pmax"],"varnames":[r"$\rm p_{max}$"],
    "ylabel":"Pa", "title":"cardiac CSF pressure"},
    {"vars":["umax"],"varnames":[r"$\rm u_{max}$"], 
    "ylabel":"mm/s", "title":"CSF velocity", "scale":1e3,},
    {"vars":["cardiac_umax"],"varnames":[r"$\rm u_{max}$"], 
    "ylabel":"m/s", "title":"cardiac CSF velocity"},
    ]

    barplot_groups += [{"vars":[f"cmean_{dom}_{t}" for t in [10800, 21600,  43200, 86400]],
    "varnames":[f"{int(t/3600)} h" for t in [10800, 21600,  43200, 86400]], 
    "ylabel":"mmol/l", "title":f"mean concentration ({dom})"} for
    dom in ["CSF", "parenchyma", "PVS artery", "PVS vein"]]

    # ==========================================================================
    # PATTERN 1: Multi-Panel Grid Spec Bar Plots (Unified CSV generation)
    # ==========================================================================
    nfigs = len(barplot_groups)
    fig_width = [ 1  + int(len(bp["vars"]) / 4) for bp in barplot_groups]
    nrows, ncols = 4,4
    fig = plt.figure(figsize=(7, 8))
    gs = plt.GridSpec(nrows, ncols, figure=fig)
    gs.update(wspace=0.9, hspace=0.9)
    sns.set_palette("crest", n_colors=len(df.columns))
    rowi, coli = 0, 0
    plot_order = range(len(barplot_groups))
    alphabet = "ABCDEFGHIJKLMNO"
    
    # Dictionary to combine all layout matrices into one shared sheet
    source_data_bars = {}

    for i in plot_order:
        bp,w = barplot_groups[i], fig_width[i]
        subdf = df.transpose()[bp["vars"]].transpose()
        subdf *= bp.get("scale", 1)
        subdf.index = bp["varnames"]
        
        # Capture panel metadata alongside configuration rows
        panel_prefix = f"Panel_{alphabet[i]}_{bp['title'].replace(' ', '_')}"
        source_data_bars[f"{panel_prefix}_metric"] = pd.Series(subdf.index.values)
        for m in df.columns:
            source_data_bars[f"{panel_prefix}_{m}"] = pd.Series(subdf[m].values)

        if coli + w > ncols: rowi += 1; coli = 0
        ax = fig.add_subplot(gs[rowi, coli:coli + w])
        bars = subdf.plot(ax=ax, kind="bar", ylabel=bp["ylabel"],
                        legend=False, title=bp["title"], rot=30)
        ax.set_ylim(bottom=0)
        ax.text(s=alphabet[i], x=-1.0, y=ax.get_ylim()[1]*1.3,fontweight="bold")
        coli += w
    
    plt.figlegend(df.columns, loc = 'outside upper center', ncol=len(models), 
                bbox_to_anchor=(0.5, 0.97), frameon=False)
    
    grid_filename = f"plots/comparisons/{v}/{v}.pdf"
    plt.savefig(grid_filename)
    plt.close()
    
    # Save Bar Charts Matrix Source Data CSV
    pd.DataFrame(source_data_bars).to_csv(grid_filename.replace(".pdf", "_source_data.csv"), index=False)
    print(f"--> Exported consolidated bar chart grid data to: {grid_filename.replace('.pdf', '_source_data.csv')}")

    # ==========================================================================
    # PATTERN 2: Boundary Undershoot Continuous Min/Max Overlays
    # ==========================================================================
    linestyles = ["dotted", "dashed", "solid"]
    for dom in ["par", "csf","art", "ven"]:
        plt.figure()
        
        # Initialize tracking dictionary for this specific domain file run
        source_data_minmax = {}

        for m,ls in zip(models, linestyles):
            tm = m
            tmconfig = read_config(f"configfiles/{tm}.yml")
            trmetrics = read_config(f"results/{tm}/{tm}_metrics.yml")
            dmax = trmetrics[f"{dom}_max"]
            dmin = trmetrics[f"{dom}_min"]
            times = np.arange(0, tmconfig["T"], tmconfig["dt"]) / (60*60)
            
            plt.plot(times, dmax, label=f"max {m}", ls=ls, color="maroon")
            plt.plot(times, dmin, label=f"min {m}", ls=ls, color="teal")
            
            # Map structural columns into safely padded rows using pd.Series
            source_data_minmax[f"{m}_time_hours"] = pd.Series(times)
            source_data_minmax[f"{m}_max_concentration"] = pd.Series(dmax)
            source_data_minmax[f"{m}_min_concentration"] = pd.Series(dmin)
            
        plt.xlabel("time (h)")
        plt.ylabel("tracer concentration (mmol/l)")
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=len(models),
                columnspacing=0.7, frameon=False)
        
        minmax_filename = f"plots/comparisons/{v}/{dom}_minmax.png"
        plt.savefig(minmax_filename)
        plt.close()
        
        # Save MinMax Source Data CSV
        pd.DataFrame(source_data_minmax).to_csv(minmax_filename.replace(".png", "_source_data.csv"), index=False)
        print(f"--> Exported undershoot line data to: {minmax_filename.replace('.png', '_source_data.csv')}")

    # ==========================================================================
    # PATTERN 3: Compartmental Tracer Concentration Means Over Time
    # ==========================================================================
    for doms, colors in zip([["par", "csf"],["art", "ven"]],
                            [["mediumvioletred","navy"], ["maroon", "teal"]]):
        plt.figure()
        marker = ["s", "o", "D"]
        
        # Initialize dictionary to track the dual domain plot outputs
        source_data_mean = {}

        for i, (m, mark, ls) in enumerate(zip(models, marker, linestyles)):
            tm = m
            tmconfig = read_config(f"configfiles/{tm}.yml")
            trmetrics = read_config(f"results/{tm}/{tm}_metrics.yml")
            times = np.arange(0, tmconfig["T"], tmconfig["dt"]) / (60*60)
            
            source_data_mean[f"{m}_time_hours"] = pd.Series(times)
            
            for dom,c in zip(doms,colors):
                dmean = trmetrics[f"{dom}_mean"]
                plt.plot(times, dmean, label=f"{dom.upper()} {m}", ls=ls, markerfacecolor="white",
                        marker=mark, color=c, ms=8, markevery=(0.05*i, 0.2))
                
                # Append individual spatial target measurements
                source_data_mean[f"{m}_{dom.upper()}_mean_concentration"] = pd.Series(dmean)
                
        plt.xlabel("time (h)")
        plt.ylabel("tracer concentration (mmol/l)")
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=len(models),
                    columnspacing=0.7, frameon=False)
        
        mean_filename = f"plots/comparisons/{v}/{v}_{'_'.join(doms)}_mean.png"
        plt.savefig(mean_filename)
        plt.close()
        
        # Save Mean Values Source Data CSV
        pd.DataFrame(source_data_mean).to_csv(mean_filename.replace(".png", "_source_data.csv"), index=False)
        print(f"--> Exported regional mean trend line data to: {mean_filename.replace('.png', '_source_data.csv')}")


if __name__ == "__main__":
    typer.run(compare_models)