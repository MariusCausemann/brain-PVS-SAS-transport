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
import joypy
from scipy.stats import binned_statistic
from cmap import Colormap
from scipy.interpolate import pchip

def plot_conc_ridgeline(model:str):

    os.makedirs(f"plots/{model}/", exist_ok=True)

    results = read_config(f"results/{model}/mean_concentrations.yml")
    config = read_config(f"configfiles/{model}.yml")
    artlabels = [l for l,p in pointlabels]
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, 3*dt*config["output_frequency"])
    set_plotting_defaults()
    labeldist = {n:results[f"{n}_root_dist"] for n in artlabels}
    xrange = (0, 0.25)
    timespoints = np.array([1,2, 3,4,5, 6, 9, 12])*3600
    timeidx = np.where(np.isin(times, timespoints))[0]
    binedges = np.linspace(*xrange, 50)
    #from IPython import embed;embed()
    pvsconc = [results["c_pvs"][ti] for ti in timeidx]
    outerconc = [results["c_averages"][ti] for ti in timeidx]
    jumpconc = [results["c_pvs"][ti] - results["c_averages"][ti] for ti in timeidx]

    pvsconc = [np.where(c > 0, c, 0) for c in pvsconc]
    outerconc = [np.where(c > 0, c, 0) for c in outerconc]
    cellvol = results["cellvol"]
    dist_to_root = results["dist_to_root"]

    binmid = 0.5*(binedges[1:] + binedges[:-1])
    xsmooth = np.linspace(*binedges[[0,-1]], 500)
    pchipsmooth = lambda x: pchip(binmid, x)(xsmooth)

    for rawdata, datatype in zip([pvsconc, outerconc, jumpconc], ["pvs", "outer", "jump"]):
        for plottype in ["total","conc"]:
            if plottype=="total":
                data, _, _ = binned_statistic(dist_to_root, [c*cellvol for c in rawdata], statistic="sum", bins=binedges)
            else:
                data, _, _ = binned_statistic(dist_to_root, rawdata, statistic="mean", bins=binedges)

            for s in ["smoothed", "raw"]:
                xdata = xsmooth if s=="smoothed" else binmid
                ridge_data = pd.DataFrame({f"{int(t/3600)} h": pchipsmooth(d) if s=="smoothed" else d 
                                        for t,d in zip(timespoints, data)}, index=xdata)
                ridge_data = ridge_data[ridge_data.columns[::-1]]
                fig, ax = joypy.joyplot(ridge_data, colormap=Colormap("Blues_r"), kind="values",
                                        x_range=(0, len(xdata)), alpha=0.8, figsize=(6,6), overlap=1.5,
                                        #title="Total PVS tracer content", 
                                        )
                ax[-1].set_xlabel("distance from root (m)")
                ax[-1].set_xticks(np.linspace(0, len(xdata), 5, endpoint=False),
                                np.linspace(*xdata[[0, -1]], 5, endpoint=False).round(2))

                secxa, secxb = ax[-1].secondary_xaxis(0.9), ax[-1].secondary_xaxis(0.9)
                toplabels, bottomlabels = artlabels[0::2], artlabels[1::2]
                for ln in ["ICA-L", "ACA-A1-L"]: toplabels.remove(ln)
                for ln in ["ACA-A1-R"]: bottomlabels.remove(ln)
                secxa.set_xticks([labeldist[ln]*len(xdata) / xdata[-1] for ln in toplabels],
                                toplabels, rotation=45,)
                secxb.set_xticks([labeldist[ln]*len(xdata) / xdata[-1] for ln in bottomlabels],
                                bottomlabels, rotation=45)
                secxb.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
                fig.savefig(f"plots/{model}/{model}_ridgeline_{datatype}_{plottype}_{s}.png",
                            bbox_inches='tight')


if __name__ == "__main__":
    typer.run(plot_conc_ridgeline)




