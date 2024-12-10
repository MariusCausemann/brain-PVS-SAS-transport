import yaml
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config
import typer
from typing import List
import os
from label_arteries import pointlabels
import pandas as pd
import seaborn as sns
import numpy as np

def plot_conc_at_label(model:str):

    os.makedirs(f"plots/{model}/", exist_ok=True)

    results = read_config(f"results/{model}/mean_concentrations.yml")
    config = read_config(f"configfiles/{model}.yml")

    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, 3*dt*config["output_frequency"])
    set_plotting_defaults()

    artlabelsets = [["BA",],#["ICA-L", "ICA-R"], 
                    ["MCA-L", "MCA-R"],["MCA2-L", "MCA2-R"], 
                    #["ACA-A1-L", "ACA-A1-R"], 
                    ["ACA-A2", "ACA-A3"],
                    #["PCA-L", "PCA-R"], 
                    #["PER-L", "PER-R"]
                    ]
    
    #from IPython import embed; embed()
    artlabels = sum(artlabelsets, [])

    coldict = dict()
    colors = ["#f72585", "#3a0ca3", "#4cc9f0"]

    for al in artlabelsets:
        coldict.update(dict(zip(al, colors)))
    print(coldict)
    pvspeaktimedict = {n:results[f"{n}_pvs_peak_time"] for n in artlabels}
    pvsftadict = {n:results[f"{n}_fta"] for n in artlabels}
    lagdict = {n:results[f"{n}_lag"] for n in artlabels}
    peak_annotate = ["BA"]
    fta_annotate = ["BA"]
    dt_annotate = ["MCA-R"]

    conc_at_point = results["conc_at_point"]
    avg_conc_around_point = results["avg_conc_around_point"]
    nrows = 1
    ncols = int(np.ceil(len(artlabelsets) / nrows))
    for annotate in [True, False]:
        fig, axes = plt.subplots(nrows, ncols, figsize=(int(ncols*3),nrows*3.5), constrained_layout=True)
        for ax, al_set in zip(axes.flatten(), artlabelsets):
            peaks, lags, ftas = [],[],[]
            for n in al_set:
                col = coldict[n]
                h1 = ax.plot(times / 3600, conc_at_point[n], color=col, label=n)
                h2 = ax.plot(times / 3600, avg_conc_around_point[n], color=col, ls="dashed")
                ax.fill_between(times / 3600, conc_at_point[n], avg_conc_around_point[n],
                                alpha=0.2, color=col)
                peaks.append(pvspeaktimedict[n]); lags.append(lagdict[n]), ftas.append(pvsftadict[n])

            if annotate:
                for n in al_set:
                    if n in peak_annotate:
                        ax.scatter(pvspeaktimedict[n] / 3600 ,  max(conc_at_point[n]), marker=10, color="black")
                        ax.annotate("peak", (pvspeaktimedict[n] / 3600 ,  max(conc_at_point[n])),
                                    horizontalalignment="left", verticalalignment="bottom",
                                    textcoords='offset pixels', xytext=(10, 0))
                    if n in fta_annotate:
                        ax.scatter(pvsftadict[n] / 3600 ,  0.1, marker="|", s=50, color="black")
                        ax.annotate("FTA", (pvsftadict[n] / 3600 ,  0.1), horizontalalignment="left",
                                    verticalalignment="top", textcoords='offset pixels', xytext=(5, 8))
                    if n in dt_annotate:    
                        ax.plot([pvspeaktimedict[n] / 3600, (pvspeaktimedict[n] + lagdict[n])/ 3600], 
                            2*[0.5*(max(conc_at_point[n]) + max(avg_conc_around_point[n]))], marker="|", color="black")
                        ax.plot([pvspeaktimedict[n] / 3600, pvspeaktimedict[n] / 3600, 
                                (pvspeaktimedict[n] + lagdict[n])/ 3600, (pvspeaktimedict[n] + lagdict[n])/ 3600 ], 
                            [max(conc_at_point[n])] 
                            +  2*[0.5*(max(conc_at_point[n]) + max(avg_conc_around_point[n]))] 
                            + [max(avg_conc_around_point[n])], "--", color="black", lw=1)

                        ax.annotate("Δt", ((pvspeaktimedict[n] + 0.5* lagdict[n])/ 3600 ,
                                0.5*(max(conc_at_point[n]) + max(avg_conc_around_point[n]))), 
                                horizontalalignment="center", verticalalignment="top",
                                textcoords='offset pixels', xytext=(-20, -10))

            format_secs = lambda s: "{:1.0f}:{:02.0f}".format(*divmod(abs(s)/60, 60))
            #format_secs = lambda s: f"{s/60} min"
            nl = chr(10)
            text = (f"peak: {('/' + nl).join(format_secs(p) for p in peaks)}h" + nl +
                    f"Δt: {('/' + nl).join(('' if l>-60 else '-') + format_secs(l) for l in lags)}h" + nl +
                    f"FTA: {('/' + nl).join(format_secs(p) for p in ftas if p is not None)}h")
            textleft =  np.mean(peaks) > 5*3600
            x_text = 0.03 if textleft else 1
            ax.set_xlim((0, 12))
            ax.text(x_text, 0.92,
                    text, transform=ax.transAxes, horizontalalignment='left' if textleft else "right",
                    verticalalignment="top")

            ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, 
                    columnspacing=0.3, frameon=False)
            ax.set_xlabel("time (h)")
            ax.set_ylabel("concentration (mmol/l)")
        leg = plt.figlegend(handles =[h1[0], h2[0]],# h3[0]],
                    labels=["PVS", "outside PVS"], 
                            #f"PVS proximity (R + {int(proximity_dist*1e3)} mm)"],
                    loc='lower center', #facecolor="white", 
                    edgecolor="black", frameon=False,
                    ncols=3, bbox_to_anchor=(0.5, 0.98))
        for lh in leg.legend_handles: lh.set_color('black')

        plt.savefig(f"plots/{model}/{model}_conc_at_label{'_annotated' if annotate else ''}.png",
                    bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    typer.run(plot_conc_at_label)




