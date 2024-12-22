import yaml
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config
import typer
from typing import List
import os
from label_arteries import pointlabels
import pandas as pd
import seaborn as sns


def plot_timings(models:List[str]):

    v = "_".join(models)
    os.makedirs(f"plots/comparisons/{v}/", exist_ok=True)

    results = dict()
    for m in models:
        print(f"reading {m}")
        md = dict()
        md.update(read_config(f"results/{m}/mean_concentrations.yml"))
        results[m] = md

    set_plotting_defaults()
    arteries = [n for n,p in pointlabels]
    ftas, outerftas, pts, otps  = [pd.DataFrame(
        {m:[results[m][f"{n}_{key}"] for n in arteries] for m in models},
                        index=arteries) for key in 
                        ["fta", "avg_fta", "pvs_peak_time", "avg_peak_time"]]

    ftas["dist"] = [results[m][f"{n}_root_dist"] for n in arteries]
    pts["dist"] = [results[m][f"{n}_root_dist"] for n in arteries]

    art_groups = {"MCA-L": ["ICA-L", "MCA-L", "MCA2-L", ], #"PER-L"],
                  "MCA-R": ["ICA-R", "MCA-R", "MCA2-R", ], # "PER-R"],
                  "PCA-L": ["BA", "PCA-L"],
                  "PCA-R": ["BA", "PCA-R"],
                  "ACA": ["BA", "ACA-A1-L", "ACA-A2", "ACA-A3",], # "ACA-A4"],
        }
    
    namedict = {"modelA": "baseline", "modelA-strongVM":"high PVS flow",
                "modelA-PVS-disp": "high PVS dispersion",
                "LargePVS":"high PVS flow (dilated)", "LargePVSA":"baseline (dilated)"}
    markerdict = {"modelA": "d", "modelA-strongVM":"s",
                "modelA-PVS-disp": "o", "LargePVS":"o", "LargePVSA":"*"}
    print(models)
    coldict = {"modelA-strongVM":"#6610f2","modelA":"#0b774d", "modelA-PVS-disp":"#f7b801",
               "LargePVS":"#f7b801", "LargePVSA":"darkred"}
    annotate_model = "LargePVSA" if "LargePVSA" in models else "modelA"
    for t, (pvsdf, outerdf) in zip(["fta", "peaktime"] , [(ftas, outerftas), (pts, otps)]):
        fig, axes = plt.subplots(ncols=len(art_groups) + 1, figsize=(9,4),
                                  width_ratios=[0]+[3]*len(art_groups), sharey=True)
        ax = axes[0]
        if t=="fta":ax.set_ylabel("first-time arrival (h)"); ax.set_ylim((0.0, 6))
        elif t=="peaktime":ax.set_ylabel("time-of-peak (h)"); ax.set_ylim((2, 9.1))
        ax.get_xaxis().set_ticks([])
        ax.spines['bottom'].set_visible(False)
        #ax.set_xlim((-0.05, 0.1))
        for i, (ax, (agn, ag))  in enumerate(zip(axes[1:], art_groups.items())):
            for m in models:
                print(m)
                group = pvsdf.loc[ag]
                outergroup = outerdf.loc[ag]
                ax.plot(group["dist"], group[m] / 3600, color=coldict[m],
                         marker=markerdict[m],
                         markersize=8, ls="dashed",
                        label=f"{namedict[m]}" if i==0 else None)
                #ax.plot(group["dist"], outergroup[m] / 3600, color=c, marker="x", markersize=5, ls="dotted",
                #        label=f"{namedict[m]} - outer PVS" if i==0 else None)
                #ax.fill_between(group["dist"], group[m] / 3600, outergroup[m] / 3600, 
                #                color=c, alpha=0.3)                     
                ax.set_title(agn, y=1.04)
                #ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17), ncol=2, 
                #        columnspacing=0.3, frameon=False)
                ax.set_xlabel("distance (m)")
                ax.spines['left'].set_visible(False)
                ax.tick_params(left=False)
                #ax.get_yaxis().set_ticks([])
                ax.set_xlim((0, 0.1))
                ax.get_xaxis().set_ticks([0,0.05, 0.1], ["0","0.05", "0.1"])
                for ln in ag:
                    ax.annotate(ln, (pvsdf["dist"].loc[ln], pvsdf[annotate_model].loc[ln] / 3600),
                                horizontalalignment="right", textcoords='offset points',
                                  xytext=(-3,8), fontsize=9)
        plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=len(models)*2,
                columnspacing=0.3, frameon=False)
        plt.tight_layout(w_pad=-0.5)
        plt.savefig(f"plots/comparisons/{v}/{v}_{t}.png", bbox_inches='tight',
                    dpi=400)

if __name__ == "__main__":
    typer.run(plot_timings)




