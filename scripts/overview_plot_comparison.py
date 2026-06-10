import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import typer
from typing import List
from plotting_utils import (time_str, read_config,
                    get_result, T1_clip_plot)
from subdomain_ids import CSFID, CSFNOFLOWID, PARID, LVID, V34ID
import pyvista as pv
from cmap import Colormap
import nibabel
from extract_vessels import np2pv
import matplotlib as mpl
from extract_vessels import get_tubes
from pyvista.core import _vtk_core as _vtk
from pyvista.core.filters import _get_output, _update_alg
from pyvista.core.utilities.helpers import generate_plane
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

ndict = {"x":(1,0,0), "y":(0,-1,0), "z":(0,0,1)}

percentile = 95

def clip_closed_surface(surf, normal='x', origin=None, tolerance=1e-06,
                        transfer_data=False):
    if isinstance(normal, str): normal = np.array(ndict[normal])
    if origin is None: origin = surf.center_of_mass()
    plane = generate_plane(normal, origin)
    collection = _vtk.vtkPlaneCollection()
    collection.AddItem(plane)
    alg = _vtk.vtkClipClosedSurface()
    alg.SetGenerateFaces(True)
    alg.SetInputDataObject(surf)
    alg.SetTolerance(tolerance)
    alg.SetClippingPlanes(collection)
    try:
        _update_alg(alg, False, 'Clipping Closed Surface')
    except TypeError:
        _update_alg(alg)
    result = _get_output(alg)
    if transfer_data: return result.sample(surf)
    return result


def subtract_fields(primary, secondary, times):
    """Replace primary's concentration with (primary - secondary), in place.

    Both datasets must come from the same mesh so that the c_{t} arrays line
    up cell-for-cell (or point-for-point). A clear error is raised otherwise.
    """
    for t in times:
        key = f"c_{t}"
        if key not in primary.array_names or key not in secondary.array_names:
            raise KeyError(f"Field '{key}' missing from one of the two models.")
        if primary[key].shape != secondary[key].shape:
            raise ValueError(
                f"Field '{key}' shape mismatch "
                f"({primary[key].shape} vs {secondary[key].shape}). "
                "The two models must share the same mesh to be differenced.")
        primary[key] -= secondary[key]
    return primary


def generate_difference_plot(model_a:str,
                             model_b:str,
                             times:List[int]):
    """Plot the difference (model_a - model_b) of the concentration fields."""
    times = np.array(times) * 3600
    cbar_title = "concentration difference (mmol/l)"

    # Geometry / structural data comes from model_a (assumed identical mesh).
    conf = read_config(f"configfiles/{model_a}.yml")
    mesh = pv.read_meshio(conf["mesh"])
    t1data = nibabel.load("data/T1.nii.gz")
    res = np.array(t1data.header["pixdim"][1:4])*1e-3
    Timg = np2pv(t1data.get_fdata()[:,:-50, :], res,
                origin=(0,0,0))

    # Load both models and difference their concentration fields.
    sas = get_result(model_a, "sas", times)
    sas_b = get_result(model_b, "sas", times)
    subtract_fields(sas, sas_b, times)

    art = get_result(model_a, "artery", times).ctp()
    art_b = get_result(model_b, "artery", times).ctp()
    subtract_fields(art, art_b, times)

    arteries = get_tubes(art)
    art["radius"] *= conf["pvs_ratio_artery"]
    art_pvs = get_tubes(art)
    print(art.array_names)

    csf = sas.extract_cells(np.isin(sas["label"],
                                    [CSFID, LVID, V34ID, CSFNOFLOWID]))
    par = sas.extract_cells(sas["label"]==PARID)
    # Linear scaling of the parenchyma field; applied to the difference
    # is equivalent to applying it to each model before subtracting.
    for t in times: par[f"c_{t}"] *= 0.2

    # Diverging colormaps and symmetric limits, since differences are signed.
    par_cmap, par_clim = "crameri:managua", (-0.2, 0.2)
    csf_cmap, csf_clim = "crameri:vanimo", (-1.0, 1.0)
    # Opacity transfer: transparent near zero (models agree), opaque at the
    # extremes (large disagreement).
    opacity = [0.85, 0.45, 0.05, 0.45, 0.85]

    styles = ['default', 'dark_background']

    for style in styles:
        plt.style.use(style)
        sns.set_context("paper", font_scale=1.5)
        fig = plt.figure(figsize=(len(times)*3, 10), frameon=True)
        grid = ImageGrid(fig, 111, nrows_ncols=(3, len(times)),
                        axes_pad=0,)
                        #cbar_location="right",cbar_mode=None,cbar_size="4%",cbar_pad=0.2)
        origin = np.array([0.0836, 0.092, 0.1])
        for j,n in enumerate(["x", "y", "z"]):
            nv = np.array(ndict[n])
            artclip = clip_closed_surface(arteries.extract_surface(), normal=-nv,
                                        origin=origin + nv*1e-4,)
            pvsclip = clip_closed_surface(art_pvs.extract_surface(), normal=-nv,
                                        origin=origin + nv*2e-4, transfer_data=True)
            T1slice = Timg.slice(normal=n, origin=origin)
            T1slice.field_data.update({"cmap":"gray"})
            parslice = par.slice(normal=n, origin=origin + nv*2e-10)
            csfslice = csf.slice(normal=n, origin=origin + nv*2e-10)
            parslice.field_data.update(dict(cmap=par_cmap, clim=par_clim, opacity=opacity))
            csfslice.field_data.update(dict(cmap=csf_cmap, clim=csf_clim, opacity=opacity))
            #csfslice.field_data.update(dict(color="green", opacity=0.5))
            #parslice.field_data.update(dict(color="red", opacity=0.5))
            pia = mesh.extract_cells(mesh["label"]==PARID).extract_surface().slice(normal=n, origin=origin + nv*2e-10)
            dura = mesh.extract_surface().slice(normal=n, origin=origin + nv*2e-10)
            pia.field_data.update(dict(color="violet",))# line_width=1.0))
            dura.field_data.update(dict(color="cyan",))# line_width=1.0))
            artclip.field_data["color"] = "red"
            pvsclip.field_data.update(csfslice.field_data)
            for i, t in enumerate(times):
                for g in [csfslice, parslice, pvsclip]: g.set_active_scalars(f"c_{t}")
                img = T1_clip_plot([T1slice, parslice, csfslice,
                                    pvsclip, artclip, pia, dura], normal=n)
                ax = grid.axes_row[j][i]
                ax.axis('off')
                ax.imshow(img)
                if j==0:
                    ax = grid.axes_row[0][i]
                    ax.text(0.5, 1.05, f"{time_str(t)} h", ha="center",
                        transform=ax.transAxes, fontsize=12)

        for name, dom, loc in zip(["CSF/PVS ", "parenchyma "], [csfslice, parslice], [1.05, 1.65]):
            fd = dom.field_data
            mplcmap = Colormap(fd["cmap"]).to_matplotlib()
            norm = mpl.colors.Normalize(*fd["clim"])
            cax = inset_axes(ax, width="10%",height="300%", loc="lower left",
            bbox_to_anchor=(loc, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0)
            plt.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=mplcmap),
                        label=name + cbar_title, extend='both', cax=cax)

        tstr = "-".join([f"{t/3600:.0f}" for t in times])
        plt.savefig(f"plots/comparisons/{model_a}_vs_{model_b}_diff_{tstr}{'_dark' if style=='dark_background' else ''}.png",
                        dpi=300, bbox_inches="tight", transparent=True)

if __name__ == "__main__":
    typer.run(generate_difference_plot)