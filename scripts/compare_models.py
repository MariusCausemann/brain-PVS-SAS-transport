import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import typer
from typing import List
from plotting_utils import (time_str, read_config, compute_ranges,
                    get_result, clip_plot, isosurf_plot, timesurf_plot,
                    compute_diff_ranges)
import pyvista as pv
from generate_synthseg_mesh import CSFID, CSFNOFLOWID, PARID, LVID, V34ID
import matplotlib as mpl

crop_dict = {"overview":lambda img : img[100:-10,140:-140],
             "overview_bottom":lambda img : img[100:-100,140:-140],
             "detail":lambda img : img,
            "isosurf":lambda img : img[100:-10,140:-140],
            "timesurf":lambda img : img[100:-10,140:-140],
            "detail_bottom":lambda img : img,
            "isosurf_bottom":lambda img : img[100:-100,140:-140],
            "timesurf_bottom":lambda img : img[100:-100,140:-140],
             }
padding = {"overview":0, "detail":0.05, "isosurf":0, "timesurf":0}
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

    model_data = dict()
    for j,m in enumerate([modela, modelb,"diff"]):
        model_data[m] = dict()
        for k in ["sas", "artery", "vein"]:
            if m=="diff":
                da = model_data[modela][k].copy()
                db = model_data[modelb][k]
                if da.n_points != db.n_points:
                    db = da.sample(db, pass_point_data=False)
                for t in times: da[f"c_{t}"] = da[f"c_{t}"] - db[f"c_{t}"]
                model_data[m][k] = da
            else:
                conf = read_config(f"configfiles/{modela}.yml")
                mesh = pv.read_meshio(conf["mesh"])
                bg = mesh.extract_surface()
                model_data[m][k] = get_result(m,k, times)
        sas = model_data[m]["sas"]
        art = model_data[m]["artery"].ctp()
        ven = model_data[m]["vein"].ctp()
        for netw in [art, ven]: netw["radius"] *= conf["pvs_ratio_artery"]
        csf = sas.extract_cells(np.isin(sas["label"], [CSFID, LVID, V34ID, CSFNOFLOWID]))
        par = sas.extract_cells(sas["label"]==PARID)
        for t in times:par[f"c_{t}"] *= 0.2
        comb = pv.merge([csf, par], merge_points=False)
        clim = (0, cmax)
        cbar_title = "concentration (mol/l)"
        for i, t in enumerate(times):
            cmap = "tempo"
            if m=="diff":
                clim=(-diffmax, diffmax)
                cbar_title = "concentration diff (mol/l)"
                cmap="curl"
            for i, t in enumerate(times):
                title = f"time: {time_str(t)} h"
                for g in [csf, art, ven, par, comb]: g.set_active_scalars(f"c_{t}")
                if type=="overview":
                    img = clip_plot(csf, par, [art, ven], None, title, clim=clim, 
                            cmap=cmap, cbar_title=cbar_title)
                if type=="isosurf":
                    img = isosurf_plot(comb, [art], bg, None, title, clim=clim, 
                            cmap=cmap, cbar_title=cbar_title)
                if type=="timesurf":
                    cmap = "colorbrewer:spectral"
                    cbar_title = f"time (h)"
                    img = timesurf_plot(comb, [art],bg, times,t, None, title, clim=clim, 
                            cmap=cmap)
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