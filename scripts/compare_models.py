import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import typer
from typing import List
from plotting_utils import time_str, read_config

def compare_models(modela:str, modelb:str, times:List[int]):

    fig = plt.figure(figsize=(5., 9.), frameon=True)
    grid = ImageGrid(fig, 111, nrows_ncols=(len(times), 3), axes_pad=0.0)
    for i, t in enumerate(times):
        for j,m in enumerate([modela, modelb, "diff"]):
            fn = f"plots/{m}/{m}_{t}.png"
            if m=="diff":
                fn = f"plots/comparisons/{modela}_{modelb}/{modela}_{modelb}_diff_{t}.png"
            img = plt.imread(fn)
            ax = grid.axes_row[i][j]
            ax.axis('off')
            if i + 1==len(times):
                ax.imshow(img[100:-10,140:-140])
            else:
                ax.imshow(img[100:-100,140:-140])

        ax = grid.axes_row[i][0]
        ax.text(-0.08, 0.5, f"{time_str(t)} h", rotation=90, va="center",
                    transform=ax.transAxes, fontsize=6)
        
    descr = [read_config(f"configfiles/{m}.yml")["description"] for m in [modela, modelb]]
        
    for j,m in enumerate(descr +  ["difference"]):
        ax = grid.axes_row[0][j]
        ax.text(0.5, 1.05, m, ha="center", transform=ax.transAxes, fontsize=6)


    plt.savefig(f"plots/comparisons/{modela}_{modelb}/{modela}_{modelb}.png",
                 dpi=300, bbox_inches="tight",)

if __name__ == "__main__":
    typer.run(compare_models)