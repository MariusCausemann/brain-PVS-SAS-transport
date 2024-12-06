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
import argparse

percentile = 95

def shrink_cbar(ax, shrink=0.9):
    b = ax.get_position()
    new_w = b.width*shrink
    pad = (b.width-new_w)/2.
    new_x0 = b.x0 + pad
    new_x1 = b.x1 - pad
    b.x0 = new_x0
    b.x1 = new_x1
    ax.set_position(b)

def get_tangent(v, n):
    return v - np.outer(np.dot(v,n), n)

def get_normal(netw, pidx):
    neighors = netw.point_neighbors(pidx)
    p1, p2 = netw.points[pidx], netw.points[neighors[0]]
    return (p1 - p2) / np.linalg.norm(p1 - p2)

def get_radius(netw, pidx):
    return netw.extract_points(pidx)["radius"].max()

def get_slice(art, sas, par, p, pvs_radius_ratio=2):
        pidx = map_kdtree(art.points, np.array([p]))[0]
        r = get_radius(art, pidx)
        n = get_normal(art, pidx)
        fp1, fp2, fp3 = p + 1e-4*n, p + 2e-4*n, p + 3e-4*n
        vessel = pv.SolidSphere(center=fp3, direction=n, outer_radius=r).slice(origin=fp3, normal=n)
        vessel.add_field_data(name="color", array="red")
        pvs = pv.SolidSphere(center=fp2, direction=n, outer_radius=r*pvs_radius_ratio).slice(origin=fp2, normal=n)
        pvs_surf = pv.Sphere(center=fp2, radius=r*pvs_radius_ratio).slice(origin=fp2, normal=n)
        pvs_surf.add_field_data(name="color", array="yellow")
        pvs_surf.add_field_data(name="line_width", array=5)
        par["n"] = get_tangent(par.point_data["Normals"], n)
        par["n"] /= np.linalg.norm(par["n"], axis=1, keepdims=True)
        pia = par.warp_by_vector("n", factor=-2e-4).slice(origin=fp1, normal=n)
        pia_inner = par.warp_by_vector("n", factor=2e-4).slice(origin=0.5*(fp1 + fp2), normal=n)
        pia_inner.add_field_data(name="color", array="salmon")
        pia.add_field_data(name="color", array="cyan")
        for k,v in art.point_data.items():
            pvs[k] = v[pidx]*np.ones(pvs.n_points)
        slice = sas.slice(origin=p, normal=n).clip_surface(pv.Sphere(center=p, radius=np.linalg.norm(4e-2*n)/3))
        return slice, pvs, vessel, pvs_surf, pia, pia_inner, n

def generate_plot(
    modelname: str,
    selected_artlabels:List[str],
    times:List[int],
    cmax: List[float],
):

    times = np.array(times)*3600
    pointdict = {l:p for l,p in pointlabels}
    if selected_artlabels is None:
        artlabels = list(pointdict.keys())
    else:
        artlabels = selected_artlabels
    points = np.array([pointdict[l] for l in artlabels])
    fontsize = 16

    fig = plt.figure(figsize=(len(times)*5, len(artlabels)*3), frameon=True)
    grid = ImageGrid(fig, 111, nrows_ncols=(len(artlabels), len(times)),
                     axes_pad=(0.02, 0.02) if cmax else (0.02, 0.22),
                     cbar_location="bottom", cbar_pad=0.03,
                     cbar_mode="edge" if cmax else "each",
                     cbar_size="5%"
                   )

    subdomains = pv.read_meshio("mesh/standard/standard.xdmf")
    par = subdomains.extract_cells(subdomains["label"]==PARID).extract_surface()
    par.compute_normals(inplace=True, cell_normals=False, consistent_normals=True,
                         non_manifold_traversal=False)
    sas = get_result(modelname, "sas", times)
    art = get_result(modelname, "artery", times)

    pvs_radius_ratio = 2
    cmap = cmocean.cm.algae_r
    for j, (p, l) in enumerate(zip(points, artlabels)):
        slice, pvs, vessel, pvs_surf, pia, pia_inner, n = get_slice(art, sas, par, p, pvs_radius_ratio=pvs_radius_ratio)
        max_val_sas = [slice[f"c_{t}"].max() for t in times]
        max_val_art = [pvs[f"c_{t}"].max() for t in times]

        for i, t in enumerate(times):
            sas.set_active_scalars(f"c_{t}")
            art.set_active_scalars(f"c_{t}")
            pvs.set_active_scalars(f"c_{t}")
            if cmax:
                vmax = cmax[i]
            else: vmax = max(max_val_sas[i], max_val_art[i], 0.1)
            cm, img = detail_plot(sas, [], [pvs, vessel, pvs_surf, pia, pia_inner], None, p,
                            clim=(0, vmax), 
                            normal=n*4e-2, cmap=cmap, cbar_title="concentration (mol/l)")
            ax = grid.axes_row[j][i]
            ax.axis('off')
            ax.imshow(img)
            norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
            if j==len(artlabels) - 1 or not cmax:
                cax = inset_axes(ax.cax, width="80%", height="80%", loc= 'lower center', bbox_to_anchor=(0, 0, 1, 1), 
                    bbox_transform=ax.cax.transAxes, borderpad=0.1)
                ax.cax.set_visible(False)
                cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                                  cax=cax, format="{x:.1f}", pad=0.03, orientation='horizontal')
                cb.ax.tick_params(labelsize=int(fontsize*0.8))
                if j==len(artlabels) -1:
                    cb.set_label("concentration (mmol/l)", fontsize=fontsize)

    for i,t in enumerate(times):
        ax = grid.axes_row[0][i]
        ax.text(0.5, 1.05, f"{time_str(t)} h", ha="center",
                    transform=ax.transAxes, fontsize=fontsize)
        
    for j,label in enumerate(artlabels):
        ax = grid.axes_row[j][0]
        ax.text( -0.13, 0.5, label, va="center",
         transform=ax.transAxes, fontsize=fontsize, rotation=90)
    tstr = "-".join([f"{t/3600:.0f}" for t in times])
    if selected_artlabels:
        artstr = "-".join(selected_artlabels) + "-"
    else: 
        artstr = ""
    if cmax:
        cmstr = "-".join([str(cm) for cm in cmax])+"_"
    else: cmstr=""

    plt.savefig(f"plots/{modelname}/{modelname}_{tstr}_{artstr}{cmstr}details.png",
                 dpi=200, transparent=True,bbox_inches="tight")


def main():
    parser = argparse.ArgumentParser(description="Process some arguments.")
    
    # Positional argument
    parser.add_argument("modelname", type=str, help="The name of the model")
    # Optional arguments
    parser.add_argument("--times", type=int, nargs='+', default=[1,3,6,12,24], help="List of times as integers")
    parser.add_argument("--cmax", type=float, nargs='+', default=None, help="List of cmax values as floats")
    parser.add_argument("--artlabels", type=str, nargs='+', default=None, help="List of art labels as strings")
    
    args = parser.parse_args()
    
    print(f"Model Name: {args.modelname}")
    print(f"Times: {args.times}")
    print(f"Cmax: {args.cmax}")
    print(f"Art Labels: {args.artlabels}")
    generate_plot(args.modelname,args.artlabels, args.times, args.cmax)

if __name__ == "__main__":
    main()