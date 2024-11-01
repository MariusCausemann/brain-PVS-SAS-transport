import pyvista as pv
import numpy as np
import os
import typer
from typing import List
from plotting_utils import get_result, time_str, clip_plot, timesurf_plot, isosurf_plot
import yaml
pv.global_theme.allow_empty_mesh = True

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

def plot_model(modelname: str, times:List[int], currenttime:int, type:str, cmax:float=None, filename:str=None):
    plotdir = f"plots/{modelname}"
    os.makedirs(plotdir, exist_ok=True)

    sas = get_result(modelname, "sas", times)
    art = get_result(modelname, "artery", times)
    ven = get_result(modelname, "vein", times)

    clim = (0, cmax) if cmax is not None else None
    if type=="overview":
        assert len(times) == 1
        t = times[0]
        title = f"time: {time_str(t)} h"
        csf = sas.extract_cells(np.isin(sas["label"], [CSFID, LVID, V34ID, CSFNOFLOWID]))
        par = sas.extract_cells(sas["label"]==PARID)
        par[f"c_{t}"] *= 0.2
        return clip_plot(csf, par, [art, ven], filename, title, clim=clim, cmap="tempo",
              cbar_title="concentration (mmol/l)")
    
    if type=="isosurf":
        assert len(times) == 1
        t = times[0]
        title = f"time: {time_str(t)} h"
        csf = sas.extract_cells(np.isin(sas["label"], [CSFID, LVID, V34ID, CSFNOFLOWID]))
        par = sas.extract_cells(sas["label"]==PARID)
        par[f"c_{t}"] *= 0.2
        csf.merge(par, inplace=True, merge_points=False)
        #sas[f"c_{t}"] *= (1 - (0.8*sas["label"]==PARID))
        return isosurf_plot(csf, [art], filename, title, clim=clim, cmap="tempo",
                            cbar_title="concentration (mmol/l)")

if __name__ == "__main__":
    typer.run(plot_model)
