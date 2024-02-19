import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import typer
from typing import List
from plotting_utils import time_str, read_config, compute_ranges, compute_diff_ranges
from generate_plot import plot_model
from create_diff_plot import plot_model_diff

crop_dict = {"overview":lambda img : img[100:-10,140:-140],
             "overview_bottom":lambda img : img[100:-100,140:-140],
             "detail":lambda img : img,
            "detail_bottom":lambda img : img,
             }
padding = {"overview":0, "detail":0.05}
percentile = 95

def compare_models(modela:str, modelb:str, type:str, cmax:float, diffmax:float, times:List[int]):

    fig = plt.figure(figsize=(5., 9.), frameon=True)
    grid = ImageGrid(fig, 111, nrows_ncols=(len(times), 3), axes_pad=padding[type])
    if not cmax:
        rangea = compute_ranges(modela, times, percentile,)
        rangeb = compute_ranges(modelb, times, percentile)
        cmax = np.max([r[1] for r in rangea.values()] + [r[1] for r in rangeb.values()])
    if not diffmax:
        rangediff = compute_diff_ranges(modela, modelb, times, percentile)
        rmin = np.min([r[0] for r in rangediff.values()])
        rmax = np.max([r[1] for r in rangediff.values()])
        diffmax = max(abs(rmin), abs(rmax))
    for i, t in enumerate(times):
        for j,m in enumerate([modela, modelb,"diff"]):
            if m=="diff":
                img = plot_model_diff(modela, modelb, t, type, cmax=diffmax)
            else:
                img = plot_model(m, t, type, cmax=cmax)
            ax = grid.axes_row[i][j]
            ax.axis('off')
            if i + 1==len(times):
                ax.imshow(crop_dict[type](img))
            else:
                ax.imshow(crop_dict[type + "_bottom"](img))

        ax = grid.axes_row[i][0]
        ax.text(-0.08, 0.5, f"{time_str(t)} h", rotation=90, va="center",
                    transform=ax.transAxes, fontsize=6)
        
    descr = [read_config(f"configfiles/{m}.yml")["description"] for m in [modela, modelb]]
        
    for j,m in enumerate(descr +  ["difference"]):
        ax = grid.axes_row[0][j]
        ax.text(0.5, 1.05, m, ha="center", transform=ax.transAxes, fontsize=6)


    plt.savefig(f"plots/comparisons/{modela}_{modelb}/{modela}_{modelb}_{type}.png",
                 dpi=300, bbox_inches="tight",)

if __name__ == "__main__":
    typer.run(compare_models)