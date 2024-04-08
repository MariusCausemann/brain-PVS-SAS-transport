import pyvista as pv
import numpy as np
import os
import typer
from typing import Tuple
from plotting_utils import get_result, time_str, clip_plot, detail_plot, isosurf_plot

def plot_model_diff(modela: str, modelb: str, t:int, type,  cmax:float=None, filename:str=None):
    plotdir = f"plots/comparisons/{modela}_{modelb}/"
    os.makedirs(plotdir, exist_ok=True)

    sasa = get_result(modela, "sas", t)
    arta = get_result(modela, "arteries", t)
    vena = get_result(modela, "venes", t)
    sasb = get_result(modelb, "sas", t)
    artb = get_result(modelb, "arteries", t)
    venb = get_result(modelb, "venes", t)
    sasa["diff"] = sasa["c_sas"] - sasb["c_sas"]
    arta["diff"] = arta["c_artery"] - artb["c_artery"]
    vena["diff"] = vena["c_vein"] - venb["c_vein"]
    sasa.set_active_scalars("diff")
    arta.set_active_scalars("diff")
    vena.set_active_scalars("diff")
    clim = (-cmax, cmax) if cmax is not None else None

    if type=="overview":
        filename = f"{plotdir}/{modela}_{modelb}_diff_{t}.png"
        title = f"time: {time_str(t)} h"
        csf = sasa.extract_cells(sasa["label"]==1)
        par = sasa.extract_cells(sasa["label"]==2)
        par["c_sas"] *= 0.2
        return clip_plot(csf, par, [arta, vena], filename, title, clim=clim, 
                  cmap="curl", cbar_title="concentration diff")
    
    if type=="detail":
        center = (0.2877, 0.17, 0.23)
        filename = f"{plotdir}/{modela}_{modelb}_diff_detail_{t}.png"
        return detail_plot(sasa, [arta, vena], filename, center, clim=clim, 
                cmap="curl", cbar_title="concentration diff")
    
    elif type=="isosurf":
        pv.global_theme.allow_empty_mesh = True
        filename = f"{plotdir}/{modela}_{modelb}_diff_isosurf_{t}.png"
        title = f"time: {time_str(t)} h"
        return isosurf_plot(sasa, [arta, vena], filename,  title, clim=clim,
                            cbar_title="concentration diff",
                            cmap="curl")

if __name__ == "__main__":
    typer.run(plot_model_diff)
