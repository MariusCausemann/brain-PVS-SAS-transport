import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import typer
from typing import List
from plotting_utils import (time_str, read_config, compute_ranges,
                    get_result, clip_plot, isosurf_plot, timesurf_plot)
import cmocean
from subdomain_ids import CSFID, CSFNOFLOWID, PARID, LVID, V34ID
import pyvista as pv
from cmap import Colormap
import nibabel
from extract_vessels import np2pv
import matplotlib as mpl

crop_dict = {"overview":lambda img : img[80:-100,200:-200],
             "overview_bottom":lambda img : img[100:-100,140:-140],
             "detail":lambda img : img,
            "isosurf":lambda img : img[100:-10,140:-140],
            "timesurf":lambda img : img[100:-10,140:-140],
            "detail_bottom":lambda img : img,
            "isosurf_bottom":lambda img : img[100:-100,140:-140],
            "timesurf_bottom":lambda img : img[100:-100,140:-140],
             }
percentile = 95

def compare_models(modela:str, modelb:str, type:str, cmax:float, times:List[int]):

    fig = plt.figure(figsize=(len(times)*3, len(times)*3), frameon=True)
    grid = ImageGrid(fig, 111, nrows_ncols=(2, len(times)),
                     cbar_location="right",
                     cbar_mode="single",
                     cbar_size="4%",
                     cbar_pad=0.1,)
    
    if not cmax:
        rangea = compute_ranges(modela, times, percentile,)
        rangeb = compute_ranges(modelb, times, percentile)
        cmax = np.max([r[1] for r in rangea.values()] + [r[1] for r in rangeb.values()])
    for j,m in enumerate([modela, modelb]):
        conf = read_config(f"configfiles/{m}.yml")
        mesh = pv.read_meshio(conf["mesh"])
        bg = mesh.extract_surface()
        t1data = nibabel.load("data/T1.nii.gz")
        res = np.array(t1data.header["pixdim"][1:4])*1e-3
        img = np2pv(t1data.get_fdata()[:,:-50, :], res, 
                    origin=(0,0,0))
        img.save("data/T1.vtk")
        bg = img.slice()
        sas = get_result(m, "sas", times)
        art = get_result(m, "artery", times).ctp()
        ven = get_result(m, "vein", times).ctp()
        for netw in [art, ven]: netw["radius"] *= conf["pvs_ratio_artery"]
        print(art.array_names)
        csf = sas.extract_cells(np.isin(sas["label"], [CSFID, LVID, V34ID, CSFNOFLOWID]))
        par = sas.extract_cells(sas["label"]==PARID)
        for t in times:par[f"c_{t}"] *= 0.2
        comb = pv.merge([csf, par], merge_points=False)
        for i, t in enumerate(times):
            title = f"time: {time_str(t)} h"
            for g in [csf, art, ven, par, comb]: g.set_active_scalars(f"c_{t}")
            cbar_title = "concentration (mol/l)"
            if type=="overview":
                cmap = "tempo_r"
                img = clip_plot(csf, par, [art, ven], None, title, clim=(0, cmax), 
                          cmap=cmap, cbar_title=cbar_title)
            if type=="isosurf":
                cmap = "tempo_r"
                img = isosurf_plot(comb, [art], bg, None, title, clim=(0, cmax), 
                          cmap=cmap, cbar_title=cbar_title)
            if type=="timesurf":
                cmap = "colorbrewer:spectral"
                cbar_title = f"time (h)"
                img = timesurf_plot(comb, [art],bg, times,t, None, title, clim=(0, cmax), 
                          cmap=cmap)
            
            ax = grid.axes_row[j][i]
            ax.axis('off')
            ax.imshow(crop_dict[type](img))

            if j==0:
                ax = grid.axes_row[0][i]
                ax.text(0.5, 1.05, f"{time_str(t)} h", ha="center",
                        transform=ax.transAxes, fontsize=12)
        
    descr = [read_config(f"configfiles/{m}.yml")["description"] for m in [modela, modelb]]
        
    for j,m in enumerate(descr):
        ax = grid.axes_row[j][0]
        ax.text( -0.08, 0.5, m, va="center", transform=ax.transAxes, fontsize=12, rotation=90)
    #from IPython import embed; embed()
    mplcmap = Colormap(cmap).to_matplotlib()
    norm = mpl.colors.Normalize(vmin=0, vmax=cmax)
    if type=="timesurf":
        norm = mpl.colors.BoundaryNorm(np.arange(0, len(times) + 1), mplcmap.N)
        cbar = grid[0].cax.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=mplcmap),
                         label=cbar_title, ticks=np.arange(0, len(times)) + 0.5)
        cbar.ax.set_yticklabels((np.array(times)/3600).astype(str))
    else:
        cbar = grid[0].cax.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=mplcmap),
                                    label=cbar_title)
    plt.savefig(f"plots/comparisons/{modela}_{modelb}/{modela}_{modelb}_{type}_horizontal.png",
                 dpi=300, bbox_inches="tight",)

if __name__ == "__main__":
    typer.run(compare_models)