import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config
import typer
from typing import List
import os
import numpy as np
from cmap import Colormap
import itertools
import seaborn as sns
import pyvista as pv
from subdomain_ids import CSFID, PARID

def plot_conc_percentiles(models:List[str]):

    group_a = ["modelA",]# "modelB1-10", "modelB1-100", "modelB1-1000"]
    group_b = ["modelA-strongVM",]# "modelB2-10", "modelB2-100", "modelB2-1000"]
    blues = sns.color_palette("Blues_d", n_colors=len(group_a))
    oranges = sns.color_palette("Oranges_d", n_colors=len(group_b))

    coldict = {n:col for n,col in zip(group_a, blues)}
    coldict.update({n:col for n,col in zip(group_b, oranges)})
    labeldict = {n:f"A - xi x {10**i}" for i,n in enumerate(group_a)}
    labeldict.update({n:f"B - xi x {10**i}" for i,n in enumerate(group_b)})
    models = group_a + group_b
    artnet = pv.read(f"results/{models[0]}/{models[0]}_artery.pvd")[1]
    marker = artnet["marker"]
    print(artnet.compute_cell_sizes().array_names)
    celllength = artnet.compute_cell_sizes()["Length"]

    q = 90
    v = "_".join(models)
    os.makedirs(f"plots/comparisons/{v}/", exist_ok=True)
    config = read_config(f"configfiles/{models[0]}.yml")
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, 3*dt*config["output_frequency"])
    results = dict()
    for m in models:
        print(f"reading model {m}...")
        md = dict()
        md.update(read_config(f"results/{m}/mean_concentrations.yml"))
        md["jump"] = [md["c_pvs"][i] - md["c_averages"][i] for i in range(len(times))]
        results[m] = md

    

    set_plotting_defaults()
    cellvol = results[models[0]]["cellvol"]
    print(f"csf cells: {(marker==CSFID).sum()}")
    print(f"par cells: {(marker==PARID).sum()}")

    cv_csf = np.where(marker==CSFID, cellvol, 0)
    cv_par = np.where(marker==PARID, cellvol, 0)
    assert np.allclose(cellvol, cv_csf + cv_par)


    plt.figure(figsize=(3,4))
    plt.bar(["CSF", "PAR"], [np.where(marker==CSFID, celllength, 0).sum(),
                             np.where(marker==PARID, celllength, 0).sum()])
    plt.ylabel("network length (m)")
    plt.xlabel("surrounding tissue type")
    plt.savefig("plots/meshplots/pvs_outer_bar.png", bbox_inches="tight")

    for weights, dom in zip([cellvol, cv_csf, cv_par], ["all", "pvs-csf", "pvs-par"]):
        for key in ["jump", "c_pvs", "c_averages"]:
            fig, ax = plt.subplots()
            #for m,col in zip(models, colors):
            #    ar = np.array(results[m][key])
                #ax.plot(times / 3600, np.percentile(ar, 50, weights=cellvol,
                #                                     axis=1, method="inverted_cdf"), color=col)
            #    ax.fill_between(times / 3600,
            #                    np.percentile(ar, q,weights=cellvol, axis=1, method="inverted_cdf"),              
            #                    np.percentile(ar, 100 - q, weights=cellvol, axis=1, method="inverted_cdf"),
            #                    alpha=0.7,color=col)
            for m in models:
                col = coldict[m]
                ar = np.array(results[m][key])
                ax.plot(times / 3600, np.average(ar, weights=weights, axis=1), color=col,
                        linewidth=2, label=labeldict[m])
            ncol=len(group_a)
            handles, labels = ax.get_legend_handles_labels()
            flip = lambda items, ncol: itertools.chain(*[items[i::ncol] for i in range(ncol)])
            ax.legend(flip(handles, ncol), flip(labels, ncol),
                    loc='upper center', bbox_to_anchor=(0.5, 1.17), ncol=ncol, 
                    columnspacing=0.3, frameon=False)        
            ax.set_xlabel("time (h)")
            ax.set_ylabel("concentration (mmol/l)")
            plt.savefig(f"plots/comparisons/{v}/{v}_percentiles_{dom}_{key}.png", 
                        bbox_inches='tight')


if __name__ == "__main__":
    typer.run(plot_conc_percentiles)




