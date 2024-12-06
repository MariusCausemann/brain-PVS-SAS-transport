import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import typer
from typing import List
from plotting_utils import (time_str, read_config,
                    get_result, T1_clip_plot)
from generate_synthseg_mesh import CSFID, CSFNOFLOWID, PARID, LVID, V34ID
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
    _update_alg(alg, False, 'Clipping Closed Surface')
    result = _get_output(alg)
    if transfer_data: return result.sample(surf)
    return result

def generate_overview_plot(model:str, 
                           times:List[int]):
    times = np.array(times) * 3600
    cbar_title = "concentration (mmol/l)"
    m = model
    conf = read_config(f"configfiles/{m}.yml")
    mesh = pv.read_meshio(conf["mesh"])
    t1data = nibabel.load("data/T1.nii.gz")
    res = np.array(t1data.header["pixdim"][1:4])*1e-3
    Timg = np2pv(t1data.get_fdata()[:,:-50, :], res, 
                origin=(0,0,0))
    sas = get_result(m, "sas", times)
    art = get_result(m, "artery", times).ctp()
    ven = get_result(m, "vein", times).ctp()
    arteries = get_tubes(art)
    for netw in [art, ven]: netw["radius"] *= conf["pvs_ratio_artery"]
    art_pvs = get_tubes(art)
    print(art.array_names)
    csf = sas.extract_cells(np.isin(sas["label"], 
                                    [CSFID, LVID, V34ID, CSFNOFLOWID]))
    par = sas.extract_cells(sas["label"]==PARID)
    for t in times:par[f"c_{t}"] *= 0.2

    styles = ['default', 'dark_background']

    for style in styles:
        plt.style.use(style)
        sns.set_context("paper", font_scale=1.5)
        fig = plt.figure(figsize=(len(times)*3, 10), frameon=True)
        grid = ImageGrid(fig, 111, nrows_ncols=(3, len(times)),
                        axes_pad=0,)
                        #cbar_location="right",cbar_mode=None,cbar_size="4%",cbar_pad=0.2)
        origin = np.array([0.0836, 0.092, 0.1])
        opacity = [0.1] + [0.8]*4
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
            parslice.field_data.update(dict(cmap="mako", clim=(0, 0.05), opacity=opacity))
            csfslice.field_data.update(dict(cmap="ember", clim=(0, 1), opacity=opacity))
            #csfslice.field_data.update(dict(color="green", opacity=0.5))
            #parslice.field_data.update(dict(color="red", opacity=0.5))
            pia = mesh.extract_cells(mesh["label"]==PARID).extract_surface().slice(normal=n, origin=origin + nv*2e-10)
            dura = mesh.extract_surface().slice(normal=n, origin=origin + nv*2e-10)
            pia.field_data.update(dict(color="violet", line_width=1))
            dura.field_data.update(dict(color="cyan", line_width=1))
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
        
        for name, dom, loc in zip(["CSF/PVS ", "parenchyma "], [csfslice, parslice], [1.05, 1.45]):
            fd = dom.field_data
            mplcmap = Colormap(fd["cmap"]).to_matplotlib()
            norm = mpl.colors.Normalize(*fd["clim"])
            #mplcmap.set_over("ivory")
            cax = inset_axes(ax, width="10%",height="300%", loc="lower left",
            bbox_to_anchor=(loc, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0)
            plt.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=mplcmap),
                        label=name + cbar_title, extend='max', cax=cax)

        tstr = "-".join([f"{t/3600:.0f}" for t in times])
        plt.savefig(f"plots/{model}/{model}_overview_{tstr}{'_dark' if style=='dark_background' else ''}.png",
                        dpi=300, bbox_inches="tight",)
    
if __name__ == "__main__":
    typer.run(generate_overview_plot)