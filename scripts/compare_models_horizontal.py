import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import typer
from typing import List
from plotting_utils import time_str, read_config, compute_ranges, compute_diff_ranges
from generate_plot import plot_model
import matplotlib as mpl
import cmocean

crop_dict = {"overview":lambda img : img[80:-100,200:-200],
             "overview_bottom":lambda img : img[100:-100,140:-140],
             "detail":lambda img : img,
            "isosurf":lambda img : img[100:-10,140:-140],
            "detail_bottom":lambda img : img,
            "isosurf_bottom":lambda img : img[100:-100,140:-140],
             }
padding = {"overview":0, "detail":0.05, "isosurf":0}
percentile = 95

def compare_models(modela:str, modelb:str, type:str, cmax:float, times:List[int]):

    fig = plt.figure(figsize=(len(times)*3, len(times)*3), frameon=True)
    grid = ImageGrid(fig, 111, nrows_ncols=(2, len(times)), axes_pad=padding[type],
                     cbar_location="right",
                     cbar_mode="single",
                     cbar_size="4%",
                     cbar_pad=0.1,)
    if not cmax:
        rangea = compute_ranges(modela, times, percentile,)
        rangeb = compute_ranges(modelb, times, percentile)
        cmax = np.max([r[1] for r in rangea.values()] + [r[1] for r in rangeb.values()])

    for i, t in enumerate(times):
        for j,m in enumerate([modela, modelb]):
            img = plot_model(m, t, type, cmax=cmax)
            ax = grid.axes_row[j][i]
            ax.axis('off')
            ax.imshow(crop_dict[type](img))

        ax = grid.axes_row[0][i]
        ax.text(0.5, 1.05, f"{time_str(t)} h", ha="center",
                    transform=ax.transAxes, fontsize=12)
        
    descr = [read_config(f"configfiles/{m}.yml")["description"] for m in [modela, modelb]]
        
    for j,m in enumerate(descr):
        ax = grid.axes_row[j][0]
        ax.text( -0.08, 0.5, m, va="center", transform=ax.transAxes, fontsize=12, rotation=90)

    cmap = cmocean.cm.tempo
    norm = mpl.colors.Normalize(vmin=0, vmax=cmax)

    grid[0].cax.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                         label="tracer concentration (mol/l)")


    plt.savefig(f"plots/comparisons/{modela}_{modelb}/{modela}_{modelb}_{type}_horizontal.png",
                 dpi=300, bbox_inches="tight",)

if __name__ == "__main__":
    typer.run(compare_models)