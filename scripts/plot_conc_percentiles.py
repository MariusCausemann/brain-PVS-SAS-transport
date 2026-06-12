import yaml
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
import pandas as pd  # <-- ADDED FOR SOURCE DATA EXPORT

flip = lambda items, ncol: itertools.chain(*[items[i::ncol] for i in range(ncol)])


def plot_conc_percentiles(models:List[str]):
    dpi = 500
    group_a = ["modelA", ]
    group_b = ["modelA-strongVM", "modelB2-10", "modelB2-100",]
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
    times = np.arange(0, T + dt, dt*config["output_frequency"])
    time_hours = times / 3600  # Baseline timeline tracking array
    
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

    # ==========================================================================
    # FIGURE Pattern 1: Network Length Bar Chart
    # ==========================================================================
    plt.figure(figsize=(3,3))
    csf_len = np.where(marker==CSFID, celllength, 0).sum()
    par_len = np.where(marker==PARID, celllength, 0).sum()
    plt.bar(["CSF", "PAR"], [csf_len, par_len])
    plt.ylabel("network length (m)")
    plt.xlabel("surrounding tissue type")
    
    bar_filename = "plots/meshplots/pvs_outer_bar.svg"
    plt.savefig(bar_filename, bbox_inches="tight", dpi=dpi)
    
    # Save Bar Chart Source Data
    pd.DataFrame({
        "Surrounding Tissue Type": ["CSF", "PAR"],
        "Network Length (m)": [csf_len, par_len]
    }).to_csv(bar_filename.replace(".svg", "_source_data.csv"), index=False)


    # ==========================================================================
    # FIGURE Pattern 2: Transient Metrics Over Time (Generates 9 loops)
    # ==========================================================================
    for weights, dom in zip([cellvol, cv_csf, cv_par], ["all", "pvs-csf", "pvs-par"]):
        for key in ["jump", "c_pvs", "c_averages"]:
            fig, ax = plt.subplots()
            
            # Initialize loop tracking dictionary
            source_data_transient = {"Time (h)": time_hours}
            
            for m in models:
                col = coldict[m]
                ar = np.array(results[m][key])
                y_avg = np.average(ar, weights=weights, axis=1)
                
                ax.plot(time_hours, y_avg, color=col, linewidth=2, label=labeldict[m])
                
                # Save out line tracks per model variation
                source_data_transient[f"{m}_mean_{key}"] = y_avg
                
            ncol = len(group_a)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(flip(handles, ncol), flip(labels, ncol),
                    loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=ncol, 
                    columnspacing=0.1, frameon=False)        
            ax.set_xlabel("time (h)")
            ax.set_ylabel("concentration (mmol/l)")
            
            img_filename = f"plots/comparisons/{v}/{v}_percentiles_{dom}_{key}.svg"
            plt.savefig(img_filename, bbox_inches='tight', dpi=dpi)
            plt.close()
            
            # Save CSV for this tracking iteration layout
            pd.DataFrame(source_data_transient).to_csv(img_filename.replace(".svg", "_source_data.csv"), index=False)
            
    # ==========================================================================
    # FIGURE Pattern 3: Parametric Metrics over Permeability
    # ==========================================================================
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
        
        # Initialize flat dict tracker for the whole parameter matrix row space
        source_data_qoi = {}

        for q, ax in zip(qois, axes):
            sns.set_palette("tab10")
            ax.set_xscale('log')
            
            for art, col in zip(["BA","MCA-R", "MCA-L", "ACA-A2", "ACA-A3"], colors):
                xis, r =  [], []
                for m in gr:
                    xis.append(results[m]["config"]["arterial_pvs_csf_permability"])
                    r.append(results[m][art + q] / (3600 if "(h)" in ylabeldict[q] else 1))
                
                ax.plot(xis, r, label=art if q==qois[0] else None, ls=":", marker="d",
                        markersize=5, color=col)
                ax.set_xlabel("PVS-CSF permeability (m/s)")
                ax.set_ylabel(ylabeldict[q])
                ax.set_title(titledict[q], y=0.98)
                
                # Format a safe label name and append coordinates
                q_clean = q.replace("_", "")
                source_data_qoi[f"{q_clean}_{art}_permeability_m_s"] = pd.Series(xis)
                source_data_qoi[f"{q_clean}_{art}_value"] = pd.Series(r)
                
        plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 1.01), ncol=len(artlabels), 
                        columnspacing=2, frameon=False)            
        
        qoi_filename = f"plots/comparisons/{va}/{va}.svg"
        plt.savefig(qoi_filename, dpi=dpi, transparent=True)
        plt.close()
        
        # Save Parametric Matrix Data
        pd.DataFrame(source_data_qoi).to_csv(qoi_filename.replace(".svg", "_source_data.csv"), index=False)


    # ==========================================================================
    # FIGURE Pattern 4: Multi-Artery Comparative Time-Series
    # ==========================================================================
    models = ["modelA-strongVM", "modelB2-100"]
    coldict = {"modelA-strongVM":"#F7B801","modelB2-100":"#5D16A6"}
    base_xi =  results["modelA-strongVM"]['config']['arterial_pvs_csf_permability']
    labeldict = {m:f"ξ x {int(results[m]['config']['arterial_pvs_csf_permability']/base_xi)}" 
                 for m in models}

    arteries = ["BA","MCA-R", "MCA-L", "ACA-A2", "ACA-A3"]
    fig, axes = plt.subplots(1, len(arteries), figsize=(3*len(arteries),3))
    plt.subplots_adjust(left=0.05, bottom=0.18, right=0.98, top=0.82, wspace=0.3, hspace=None)

    # Initialize tracker with shared timeframe baseline column
    source_data_multi_art = {"Time (h)": time_hours}

    for art, ax in zip(arteries, axes):
        for m in models:
            pvs_vals = results[m]["conc_at_point"][art]
            outer_vals = results[m]["avg_conc_around_point"][art]
            
            ax.plot(time_hours, pvs_vals, color=coldict[m], lw=2, ls="dashed",
                    label=labeldict[m] + " (PVS)" if art==arteries[0] else None)
            ax.plot(time_hours, outer_vals, ls=":", color=coldict[m], lw=2,
                    label=labeldict[m] + " (outer)" if art==arteries[0] else None)
            
            ax.fill_between(time_hours, outer_vals, pvs_vals, color=coldict[m], alpha=0.3)
            
            # Record structural line trajectory arrays
            source_data_multi_art[f"{art}_{m}_PVS_concentration_mmol_l"] = pvs_vals
            source_data_multi_art[f"{art}_{m}_outer_concentration_mmol_l"] = outer_vals

        ncol = len(group_a)
        handles, labels = ax.get_legend_handles_labels()
        plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 1.01), ncol=4, 
                columnspacing=2, frameon=False)        
        ax.set_xlabel("time (h)")
        ax.set_ylabel("concentration (mmol/l)")
        ax.set_title(art, y=0.98)
        ax.set_xlim(0, 12)

    multi_art_filename = f"plots/comparisons/{v}/{v}_conc.svg"
    plt.savefig(multi_art_filename, transparent=True, dpi=dpi)
    plt.close()
    
    # Save Final Comparative Time-Series Data Frame
    pd.DataFrame(source_data_multi_art).to_csv(multi_art_filename.replace(".svg", "_source_data.csv"), index=False)



if __name__ == "__main__":
    typer.run(plot_conc_percentiles)