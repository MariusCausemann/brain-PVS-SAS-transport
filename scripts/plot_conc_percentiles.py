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
from label_arteries import pointlabels
from subdomain_ids import CSFID, PARID

flip = lambda items, ncol: itertools.chain(*[items[i::ncol] for i in range(ncol)])


def plot_conc_percentiles(models:List[str]):
    dpi = 500
    group_a = ["modelA", ]#"modelB1-10", "modelB1-100", "modelB1-1000"]
    group_b = ["modelA-strongVM", "modelB2-10", "modelB2-100",]# "modelB2-1000"]
    blues = sns.color_palette("Blues_d", n_colors=len(group_a))
    oranges = sns.color_palette("Oranges_d", n_colors=len(group_b))
    artlabels = [l for l,p in pointlabels]
    coldict = {n:col for n,col in zip(group_a, blues)}
    coldict.update({n:col for n,col in zip(group_b, oranges)})
    labeldict = {n:f"BL - ξx{10**i}" for i,n in enumerate(group_a)}
    labeldict.update({n:f"BL+VM - ξx{10**i}" for i,n in enumerate(group_b)})
    models = group_a + group_b
    artnet = pv.read(f"results/{models[0]}/{models[0]}_artery.pvd")[1]
    marker = artnet["marker"]
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
        md["config"] = read_config(f"configfiles/{m}.yml")
        for n in artlabels:
            md[f"{n}_c_peak"] = max(md["conc_at_point"][n])
            md[f"{n}_delta_c_peak"] = max(md["conc_at_point"][n]) - max(md["avg_conc_around_point"][n])
        md["jump_at_point"] = {n:np.array(md["conc_at_point"][n]) 
                               - np.array(md["avg_conc_around_point"][n]) for n in artlabels}

        results[m] = md

    set_plotting_defaults()
    cellvol = results[models[0]]["cellvol"]
    cv_csf = np.where(marker==CSFID, cellvol, 0)
    cv_par = np.where(marker==PARID, cellvol, 0)
    assert np.allclose(cellvol, cv_csf + cv_par)

    # plot distribution of surrounding tissue
    plt.figure(figsize=(3,3))
    plt.bar(["CSF", "PAR"], [np.where(marker==CSFID, celllength, 0).sum(),
                             np.where(marker==PARID, celllength, 0).sum()])
    plt.ylabel("network length (m)")
    plt.xlabel("surrounding tissue type")
    plt.savefig("plots/meshplots/pvs_outer_bar.png", bbox_inches="tight",dpi=dpi)


    # plot mean PVS, outer PVS, and jump over time
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
            ax.legend(flip(handles, ncol), flip(labels, ncol),
                    loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=ncol, 
                    columnspacing=0.1, frameon=False)        
            ax.set_xlabel("time (h)")
            ax.set_ylabel("concentration (mmol/l)")
            plt.savefig(f"plots/comparisons/{v}/{v}_percentiles_{dom}_{key}.png", 
                        bbox_inches='tight',dpi=dpi)
            plt.close()
            
    qois = [f"_fta", f"_lag", f"_pvs_peak_time", f"_c_peak", f"_delta_c_peak"]
    ylabeldict = {f"_fta" : f"FTA (h)", f"_lag" : f"Δt (h)",
                            f"_pvs_peak_time": f"peak time (h)", 
                            f"_c_peak": r"$c_{\rm peak}$ (mmol/l)",
                            f"_delta_c_peak":r"$\Delta c_{\rm peak}$  (mmol/l)"}
    titledict = {f"_fta" : f"first-time arrival", f"_lag" : f"peak time difference",
                            f"_pvs_peak_time": f"PVS peak time", 
                            f"_c_peak": "peak PVS concentration",
                            f"_delta_c_peak":"peak concentration difference"}
    colors = ['#5d16a6','#9683ec','#f7b801','#f18701','#f35b04']
    for gr in [group_b]:
        va = "_".join(gr)
        os.makedirs(f"plots/comparisons/{va}/", exist_ok=True)
        fig, axes = plt.subplots(1, len(qois), figsize=(3*len(qois),3))
        plt.subplots_adjust(left=0.05, bottom=0.18, right=0.98, top=0.82, wspace=0.3, hspace=None)
        for q,ax in zip(qois, axes):
            sns.set_palette("tab10")
            ax.set_xscale('log')
            # plot qois over xi
            for art, col in zip(["BA","MCA-R", "MCA-L", "ACA-A2", "ACA-A3"], colors):
                xis, r =  [], []
                for m in gr:
                    xis.append(results[m]["config"]["arterial_pvs_csf_permability"])
                    r.append(results[m][art + q] / (3600 if "(h)" in ylabeldict[q] else 1))
                ax.plot(xis,r, label=art if q==qois[0] else None, ls=":", marker="d",
                        markersize=5, color=col)
                ax.set_xlabel("PVS-CSF permeability (m/s)")
                ax.set_ylabel(ylabeldict[q])
                ax.set_title(titledict[q], y=0.98)
        plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 1.01), ncol=len(art), 
                        columnspacing=2, frameon=False)            
        plt.savefig(f"plots/comparisons/{va}/{va}.png", dpi=dpi, transparent=True,
                )
        plt.close()


    #from IPython import embed; embed()
    models = ["modelA-strongVM", "modelB2-100"]
    coldict = {"modelA-strongVM":"#F7B801","modelB2-100":"#5D16A6"}
    conc_types = {"conc_at_point":"pvs", "avg_conc_around_point":"outer",
                  "jump_at_point":"jump"}   
    base_xi =  results["modelA-strongVM"]['config']['arterial_pvs_csf_permability']
    labeldict = {m:f"ξ x {int(results[m]['config']['arterial_pvs_csf_permability']/base_xi)}" 
                 for m in  models}

    arteries = ["BA","MCA-R", "MCA-L", "ACA-A2", "ACA-A3"]
    fig, axes = plt.subplots(1, len(arteries), figsize=(3*len(arteries),3))
    plt.subplots_adjust(left=0.05, bottom=0.18, right=0.98, top=0.82, wspace=0.3, hspace=None)

    for art, ax in zip(arteries, axes):
        for m in models:
            ax.plot(times / 3600, results[m]["conc_at_point"][art],
                    color=coldict[m], lw=2, ls="dashed",
                    label=labeldict[m] + " (PVS)" if art==arteries[0] else None)
            ax.plot(times / 3600, results[m]["avg_conc_around_point"][art],ls=":", 
                    color=coldict[m], lw=2,
                      label=labeldict[m] + " (outer)" if art==arteries[0] else None)
            #ax.plot(times / 3600, results[m]["avg_conc_nearby_point"][art],ls="-.", 
            #        color=coldict[m], lw=2, label=labeldict[m] + " (nearby)")
            ax.fill_between(times / 3600, results[m]["avg_conc_around_point"][art],
                             results[m]["conc_at_point"][art],
                             color=coldict[m], alpha=0.3)
        ncol=len(group_a)
        handles, labels = ax.get_legend_handles_labels()
        plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 1.01), ncol=4, 
                columnspacing=2, frameon=False)        
        ax.set_xlabel("time (h)")
        ax.set_ylabel("concentration (mmol/l)")
        ax.set_title(art, y=0.98)
        ax.set_xlim(0, 12)

    plt.savefig(f"plots/comparisons/{v}/{v}_conc.png", 
                transparent=True,dpi=dpi)
    plt.close()

    from IPython import embed; embed()



if __name__ == "__main__":
    typer.run(plot_conc_percentiles)



