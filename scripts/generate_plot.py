import pyvista as pv
import numpy as np
import os
import typer
from typing import Tuple
from plotting_utils import get_result, time_str, clip_plot, detail_plot
import yaml


def plot_model(modelname: str, t:int, type, cmax:float=None, filename:str=None):
    plotdir = f"plots/{modelname}"
    os.makedirs(plotdir, exist_ok=True)

    sas = get_result(modelname, "sas", t)
    art = get_result(modelname, "arteries", t)
    ven = get_result(modelname, "venes", t)
    clim = (0, cmax) if cmax is not None else None
    if type=="overview":
        title = f"time: {time_str(t)} h"
        return clip_plot(sas, [art, ven], filename, title, clim=clim, cmap="matter",
              cbar_title="concentration")
    
    elif type=="detail":
        center = (0.2877, 0.17, 0.23)
        return detail_plot(sas, [art, ven], filename, center, clim=clim, cmap="matter",
                cbar_title="concentration")

if __name__ == "__main__":
    typer.run(plot_model)
