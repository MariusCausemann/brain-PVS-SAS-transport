from plotting_utils import get_result, read_config
import pyvista as pv
import typer
import numpy as np
from pathlib import Path
from tqdm import tqdm
from extract_vessels import get_tubes
import k3d.colormaps.paraview_color_maps as pcm
import k3d
import matplotlib
import skimage

hexcolor = lambda c: int(matplotlib.colors.to_hex(c)[1:], base=16)
compression_level = 6

def make_html(modelname:str, n:int):
    config = read_config(f"configfiles/{modelname}.yml")
    clim = (0,2)
    #times = (np.array([0,1, 3, 6, 12])*3600).astype(np.int32)
    times = (np.array([0, 1/3, 2/3, 1,2, 3,4,5,6,8,10,12,15,18,24])*3600).astype(np.int32)
    meshfile = Path(config["mesh"])
    skull = pv.read(meshfile.parent / "surfaces/skull.ply")
    LV = pv.read(meshfile.parent / "surfaces/LV.ply")
    V34 = pv.read(meshfile.parent / "surfaces/V34.ply")
    parenchyma = pv.read(meshfile.parent / "surfaces/parenchyma_incl_ventr.ply")
    art = get_result(modelname, "artery", times).ctp()
    sas = get_result(modelname, "sas", times)
    arteries = get_tubes(art)
    arteries["c"] = arteries["c_0"]
    for t in times:
        sas[f"c_{t}"] = sas[f"c_{t}"].clip(min=0).astype(np.float32)
            
    img = pv.create_grid(sas, dimensions=(n, n, n))
    img = img.sample(sas, progress_bar=True)
    for t in times:
        img[f"c_{t}"] = skimage.filters.gaussian(img[f"c_{t}"], sigma=1, truncate=2)

    pl = k3d.plot(
        camera_rotate_speed=3,
        camera_zoom_speed=5,
        background_color=hexcolor("black"),
        grid_visible=False,
        camera_auto_fit=True,
        axes_helper=False,
        lighting=3,
        )

    vol_dtype = np.float16
    cmap = np.array(pcm.Yellow___Gray___Blue, dtype=vol_dtype).reshape(-1,4)[:,:] 
    cmap[:,0] = np.linspace(-1, 1, cmap.shape[0])[::-1]
    cmap = cmap.flatten()
    get_c = lambda t: img[f"c_{t}"].clip(min=0.0).astype(vol_dtype).reshape([n]*3)
    get_art_c = lambda t: arteries[f"c_{t}"].clip(min=0.0).astype(vol_dtype)

    pl_volume = k3d.volume({str(t/(60*60)): get_c(t) for t in times},
                           bounds=img.bounds, color_range=clim, color_map=cmap,
                           compression_level=compression_level,
                           opacity_function=np.array([(0, 0.2),(0.001, 0.8), (0.5, 0.8), (1, 0.7)], dtype=vol_dtype), 
                           alpha_coef=50, gradient_step=0.1,
                           name="concentration field")
    pl += pl_volume
    for surf, name in zip([skull, LV, V34, parenchyma],
                           ["arachnoid membrane", "lateral ventricle",
                             "3rd & 4th ventricle","parenychyma"]):
        pl_surf = k3d.vtk_poly_data(surf, wireframe=True, name=name, 
                                 color=hexcolor("silver"), opacity=0.05, 
                                 compression_level=compression_level)
        if name=="parenychyma": pl_surf.visible = False
        pl += pl_surf

    arteries = arteries.extract_geometry()
    pl_art = k3d.vtk_poly_data(arteries, color_range=clim, name="arterial PVS",
                                color_map=cmap, compression_level=compression_level)
    pl_art.attribute = {str(t/(60*60)): get_art_c(t) for t in times}
    pl += pl_art
    #sas = sas.extract_geometry()
    #pl_sas = k3d.vtk_poly_data(sas, color_range=clim, name="CSF-PAR", color_map=cmap,
    #                            compression_level=compression_level)
    #pl_sas.attribute = {str(t/(60*60)): get_sas_c(t) for t in times}
    #pl += pl_sas

    timetxt = k3d.text2d(f"", size=4, label_box=False, name="time",
                          color=hexcolor("white"), is_html=True)
    timetxt.text = {str(t):f"time: {t} h" for t in range(int(times[-1] / 3600))}
    pl += timetxt

    with open(f"plots/{modelname}/{modelname}_{n}.html", 'w') as f:
        f.write(pl.get_snapshot())


if __name__ == "__main__":
    typer.run(make_html)

