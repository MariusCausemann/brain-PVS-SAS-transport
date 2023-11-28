import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import typer
from typing import List

def compare_models(modela:str, modelb:str, times:List[int]):

    fig = plt.figure(figsize=(5., 9.), frameon=True)
    grid = ImageGrid(fig, 111, nrows_ncols=(len(times), 2), axes_pad=0.0)
    for i, t in enumerate(times):
        for j,m in enumerate([modela, modelb]):
            fn = f"plots/{m}/{m}_{t}.png"
            img = plt.imread(fn)
            ax = grid.axes_column[j][i]
            ax.axis('off')
            ax.imshow(img)

    plt.savefig(f"plots/comparisons/{modela}_{modelb}/{modela}_{modelb}.png",
                 dpi=300, bbox_inches="tight",)

if __name__ == "__main__":
    typer.run(compare_models)