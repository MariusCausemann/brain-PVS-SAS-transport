import pyvista as pv
import numpy as np
import os
import typer
from plotting_utils import get_result, time_str, clip_plot

def plot_model(modelname: str, t:int):
    plotdir = f"plots/{modelname}"
    os.makedirs(plotdir, exist_ok=True)

    sas = get_result(modelname, "sas", t)
    art = get_result(modelname, "arteries", t)
    ven = get_result(modelname, "venes", t)
    filename = f"{plotdir}/{modelname}_{t}.png"
    title = f"time: {time_str(t)} h"
    clip_plot(sas, [art, ven], filename, title, clim=(0, 1), cmap="matter",
              cbar_title="concentration")

if __name__ == "__main__":
    typer.run(plot_model)
