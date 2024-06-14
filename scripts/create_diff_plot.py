import pyvista as pv
import numpy as np
import os
import typer
from typing import Tuple
from plotting_utils import get_result, time_str, clip_plot, detail_plot, isosurf_plot

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

def plot_model_diff(modela: str, modelb: str, t:int, type,  cmax:float=None, filename:str=None):
    plotdir = f"plots/comparisons/{modela}_{modelb}/"
    os.makedirs(plotdir, exist_ok=True)

    sasa = get_result(modela, "sas", t)
    arta = get_result(modela, "artery", t)
    vena = get_result(modela, "vein", t)
    sasb = get_result(modelb, "sas", t)
    artb = get_result(modelb, "artery", t)
    venb = get_result(modelb, "vein", t)
    sasa["diff"] = sasa[f"c_{t}"] - sasb[f"c_{t}"]
    arta["diff"] = arta[f"c_{t}"] - artb[f"c_{t}"]
    vena["diff"] = vena[f"c_{t}"] - venb[f"c_{t}"]
    sasa.set_active_scalars("diff")
    arta.set_active_scalars("diff")
    vena.set_active_scalars("diff")
    clim = (-cmax, cmax) if cmax is not None else None

    if type=="overview":
        filename = f"{plotdir}/{modela}_{modelb}_diff_{t}.png"
        title = f"time: {time_str(t)} h"
        csf = sasa.extract_cells(np.isin(sasa["label"], [CSFID, LVID, V34ID]))
        par = sasa.extract_cells(sasa["label"]==PARID)
        par[f"c_{t}"] *= 0.2
        return clip_plot(csf, par, [arta, vena], filename, title, clim=clim, 
                  cmap="curl", cbar_title="concentration diff (mmol/l)")
    
    if type=="detail":
        center = (0.2877, 0.17, 0.23)
        filename = f"{plotdir}/{modela}_{modelb}_diff_detail_{t}.png"
        return detail_plot(sasa, [arta, vena], filename, center, clim=clim, 
                cmap="curl", cbar_title="concentration diff (mmol/l)")
    
    elif type=="isosurf":
        pv.global_theme.allow_empty_mesh = True
        filename = f"{plotdir}/{modela}_{modelb}_diff_isosurf_{t}.png"
        title = f"time: {time_str(t)} h"
        return isosurf_plot(sasa, [arta, vena], filename,  title, clim=clim,
                            cbar_title="concentration diff (mmol/l)",
                            cmap="curl")

if __name__ == "__main__":
    typer.run(plot_model_diff)
