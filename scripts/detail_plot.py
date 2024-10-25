import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import typer
from typing import List
from plotting_utils import time_str, read_config, detail_plot, get_result
from generate_plot import plot_model
import matplotlib as mpl
import cmocean
from label_arteries import pointlabels
from test_map_on_global_coords_shift import map_kdtree
import pyvista as pv
from generate_synthseg_mesh import PARID
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

percentile = 95

def get_tangent(v, n):
    return v - np.outer(np.dot(v,n), n)

def get_normal(netw, pidx):
    neighors = netw.point_neighbors(pidx)
    p1, p2 = netw.points[pidx], netw.points[neighors[0]]
    return (p1 - p2) / np.linalg.norm(p1 - p2)

def get_radius(netw, pidx):
    return netw.extract_points(pidx)["radius"].max()

def compare_models(modelname:str):
    times = np.array([1, 3, 6, 12, 24])*3600
    n_arteries = len(pointlabels)
    points = np.array([c for n,c in pointlabels[:n_arteries]])
    artlabels = [n for n,c in pointlabels[:n_arteries]]
    fig = plt.figure(figsize=(len(artlabels)*4, len(times)*3), frameon=True)
    grid = ImageGrid(fig, 111, nrows_ncols=(len(artlabels), len(times)),
                     axes_pad=(0.27, 0.03),
                    #cbar_location="right",
                    #cbar_mode="each",
                   )

    subdomains = pv.read_meshio("mesh/standard/standard.xdmf")
    par = subdomains.extract_cells(subdomains["label"]==PARID).extract_surface()
    par.compute_normals(inplace=True, cell_normals=False, consistent_normals=True,
                         non_manifold_traversal=False)
    sas = get_result(modelname, "sas", times)
    art = get_result(modelname, "artery", times)

    idx = map_kdtree(art.points, points)
    normals = [get_normal(art, pidx) for pidx in idx]
    radii = [get_radius(art, pidx) for pidx in idx]

    pvs_radius_ratio = 2

    cmap = cmocean.cm.algae_r
    for j, (pidx, p, l, n, r) in enumerate(zip(idx, points, artlabels, normals, radii)):
        fp1, fp2 = p + 1e-4*n, p + 2e-4*n
        vessel = pv.SolidSphere(center=fp2, direction=n, outer_radius=r).slice(origin=fp2, normal=n)
        vessel.add_field_data(name="color", array="red")
        pvs = pv.SolidSphere(center=fp1, direction=n, outer_radius=r*pvs_radius_ratio).slice(origin=fp1, normal=n)
        pvs_surf = pv.Sphere(center=fp2, radius=r*pvs_radius_ratio).slice(origin=fp2, normal=n)
        pvs_surf.add_field_data(name="color", array="yellow")
        pvs_surf.add_field_data(name="line_width", array=5)
        par["n"] = get_tangent(par.point_data["Normals"], n)
        par["n"] /= np.linalg.norm(par["n"], axis=1, keepdims=True)
        pia = par.warp_by_vector("n", factor=2e-4).slice(origin=fp1, normal=n)
        pia_inner = par.warp_by_vector("n", factor=-2e-4).slice(origin=fp2, normal=n)
        pia_inner.add_field_data(name="color", array="salmon")
        pia.add_field_data(name="color", array="cyan")

        slice = sas.slice(origin=p, normal=n).clip_surface(pv.Sphere(center=p, radius=np.linalg.norm(4e-2*n)/3))
        max_val_sas = [slice[f"c_{t}"].max() for t in times]
        max_val_art = [art[f"c_{t}"][pidx] for t in times]

        for i, t in enumerate(times):
            pvs[f"c_{t}"] = art[f"c_{t}"][pidx]*np.ones(pvs.n_points)
            sas.set_active_scalars(f"c_{t}")
            art.set_active_scalars(f"c_{t}")
            pvs.set_active_scalars(f"c_{t}")
            vmax = max(max_val_sas[i], max_val_art[i], 0.01)
            cm, img = detail_plot(sas, [], [pvs, vessel, pvs_surf, pia, pia_inner], None, p,
                            clim=(0, vmax), 
                            normal=n*4e-2, cmap=cmap, cbar_title="concentration (mol/l)")
            ax = grid.axes_row[j][i]
            ax.axis('off')
            ax.imshow(img)
            norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
            cax = inset_axes(ax, width="3%", height="60%", loc=3, bbox_to_anchor=(1.09,0.2, 1,1), 
                 bbox_transform=ax.transAxes, borderpad=0.1)

            cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),cax=cax, format="{x:.1f}",
                             shrink=0.6, pad=0.03)
            bottomtick = "0\nmol/l" if i==len(times) -1 else "0"
            cb.ax.text(0.5, -0.02, bottomtick, fontsize=8, transform=cb.ax.transAxes, va='top', ha='center')
            cb.ax.text(0.5, 1.0, f'{vmax:.6f}'[:4], fontsize=8, transform=cb.ax.transAxes, va='bottom', ha='center')
            cb.set_ticks([])
            #cb.ax.tick_params(labelsize=8)

    for i,t in enumerate(times):
        ax = grid.axes_row[0][i]
        ax.text(0.5, 1.05, f"{time_str(t)} h", ha="center",
                    transform=ax.transAxes, fontsize=12)
        
    for j,label in enumerate(artlabels):
        ax = grid.axes_row[j][0]
        ax.text( -0.13, 0.5, label, va="center",
         transform=ax.transAxes, fontsize=12, rotation=90)

    #grid[0].cax.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    #                     label="tracer concentration (mol/l)")
    #fig.tight_layout(h_pad=1.0, w_pad=1.2)

    plt.savefig(f"plots/{modelname}/{modelname}_details.png",
                 dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    typer.run(compare_models)