import pyvista as pv
import numpy as np
import os
import typer
from plotting_utils import get_result, time_str, clip_plot

def plot_model_diff(modela: str, modelb: str, t:int):
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

    filename = f"{plotdir}/{modela}_{modelb}_diff_{t}.png"
    title = f"time: {time_str(t)} h"
    clip_plot(sasa, [arta, vena], filename, title, clim=(-0.5, 0.5), 
              cmap="curl", cbar_title="concentration diff")

if __name__ == "__main__":
    typer.run(plot_model_diff)
