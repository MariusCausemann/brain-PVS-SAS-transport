from plotting_utils import get_result, read_config
import pyvista as pv
import typer
import numpy as np
from pathlib import Path
from tqdm import tqdm
from plot_csf_flow import from_k3d
import k3d.colormaps.paraview_color_maps as pcm
from extract_vessels import get_tubes
import seaborn as sns
import os
from cmap import Colormap

def make_3D_time_view(modelname:str):
    config = read_config(f"configfiles/{modelname}.yml")
    os.makedirs(f"plots/{modelname}/3D_FTA/", exist_ok=True)
    nres = 300
    times = np.array([1, 3, 6, 9, 12]) * 3600
    clim = (0, times[-1] / 3600)
    scalarbar_args = dict(color="black", title="First-time arrival (h)", 
                          vertical=False, height=0.08, width=0.6,
                          position_x=0.2, position_y=0.05, 
                          title_font_size=36, use_opacity=False,
                          label_font_size=28, fmt="%.1f")
    meshfile = Path(config["mesh"])
    skull = pv.read(meshfile.parent / "surfaces/skull.ply")
    ventricles = pv.read(meshfile.parent / "surfaces/ventricles.ply")
    art = get_result(modelname, "artery", times).ctp()
    ven = get_result(modelname, "vein", times).ctp()
    sas = get_result(modelname, "sas", times)
    arteries = get_tubes(art)
    veins = get_tubes(ven)

    img = pv.ImageData(pv.create_grid(sas, dimensions=[nres]*3))
    img = img.sample(sas, progress_bar=True)
    img["FTA"] = np.full(img[f"c_{times[0]}"].shape, times[-1] + 10)
    arteries["FTA"] = np.full(arteries[f"c_{times[0]}"].shape, np.nan)
    veins["FTA"] = np.full(veins[f"c_{times[0]}"].shape, np.nan)
    for t in times[::-1]:
        print(t)
        img["FTA"] = np.where(img[f"c_{t}"] > 0.1, t/3600, img["FTA"])
        arteries["FTA"] = np.where(arteries[f"c_{t}"] > 0.1, t/3600, arteries["FTA"])
        veins["FTA"] = np.where(veins[f"c_{t}"] > 0.1, t/3600, veins["FTA"])

    pl = pv.Plotter(off_screen=True, window_size=(1200, 1200))
    pl.set_background("white")

    cmap = from_k3d(pcm.Blue___Green___Orange).reversed()
    #cmap = Colormap("magma_r").to_matplotlib()

    pl.add_mesh(arteries, scalars="FTA", cmap=cmap,
                clim=clim, show_scalar_bar=False,
                nan_opacity=0.3, nan_color="silver")
    
    pl.add_mesh(veins, scalars="FTA", cmap=cmap,
                clim=clim, show_scalar_bar=False,
                nan_opacity=0.3, nan_color="silver")

    pl.add_volume(img,
                scalars="FTA", specular=0.8, diffuse=1, 
                ambient=0.5, blending="composite",
                cmap=cmap, opacity=[0.35]*100 + [0], clim=clim,
                show_scalar_bar=True,
                scalar_bar_args=scalarbar_args,
                opacity_unit_distance=0.002,)
    
    pl.camera_position = 'xz'
    pl.camera.azimuth += 140
    pl.camera.zoom(1.3)
    pl.screenshot(f"plots/{modelname}/3D_FTA/3D_FTA.png", 
                    transparent_background=False)
    pl.close()

    for netw, netw_name in zip([arteries, veins],["arteries", "veins"]): 
        pl = pv.Plotter(off_screen=True, window_size=(1200, 1200))
        pl.set_background("white")
        pl.add_mesh(netw, scalars="FTA", cmap=cmap,
                    clim=clim, show_scalar_bar=False,
                    nan_opacity=0.3, nan_color="silver")
        
        pl.camera_position = 'xz'
        pl.camera.azimuth += 140
        pl.camera.zoom(1.3)
        pl.screenshot(f"plots/{modelname}/3D_FTA/3D_FTA_{netw_name}.png", 
                        transparent_background=False)
        pl.close()


    
    img.save("FTA.vtk")


if __name__ == "__main__":
    typer.run(make_3D_time_view)

