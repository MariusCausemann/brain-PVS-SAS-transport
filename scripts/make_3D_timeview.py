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

def make_3D_time_view(modelname:str):
    config = read_config(f"configfiles/{modelname}.yml")
    os.makedirs(f"plots/{modelname}/timeview3D/", exist_ok=True)
    colors = [ "#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"]
    getcmap = lambda color: sns.blend_palette(["white", color],as_cmap=True)
    clim = (0,2)
    nres = 350
    times = np.array([1, 6, 12, 24]) * 3600
    
    meshfile = Path(config["mesh"])
    skull = pv.read(meshfile.parent / "surfaces/skull.ply")
    ventricles = pv.read(meshfile.parent / "surfaces/ventricles.ply")
    art = get_result(modelname, "artery", times).ctp()
    sas = get_result(modelname, "sas", times)
    arteries = get_tubes(art)
    for t in times:
        sas[f"c_{t}"] = np.where(sas[f"c_{t}"] > clim[1], 
                                 clim[1], sas[f"c_{t}"])
    img = pv.ImageData(pv.create_grid(sas, dimensions=[nres]*3))
    img = img.sample(sas, progress_bar=True)

    for type in ["volume", "contour"]:
        pl = pv.Plotter(off_screen=True, window_size=(1200, 1200))
        pl.set_background("white")
        pl.add_mesh(skull,style='points',point_size=0.1, opacity=0.6, show_scalar_bar=False,
                    color="silver", specular=1, render_points_as_spheres=True)
        pl.add_mesh(ventricles, style='points',point_size=0.8,
                    show_scalar_bar=False,specular=1,
                    render_points_as_spheres=True,
                    color="silver")
        #pl.add_mesh(arteries, color="silver",
        #            clim=clim, show_scalar_bar=False,)
        pl.camera_position = 'xz'
        pl.camera.azimuth += 140
        pl.camera.zoom(1.5)

        # generate frames
        for t, color in zip(times, colors):
            cmap = getcmap(color)
            pl.add_mesh(arteries, scalars=f"c_{t}", cmap=cmap,
                        clim=clim, show_scalar_bar=False,)
            if type=="volume":
                pl.add_volume(img,
                          scalars=f"c_{t}", specular=0.8, diffuse=1, 
                          ambient=0.5, blending="composite",
                          cmap=cmap, opacity=[0, 0.5], clim=clim,
                          show_scalar_bar=False,
                          opacity_unit_distance=0.005,)
            elif type=="contour":
                pl.add_mesh(sas.contour(np.linspace(clim[-1]/2, clim[-1], 5),
                            scalars=f"c_{t}").smooth(100),
                            scalars=f"c_{t}", opacity=0.5,
                            clim=clim, cmap=cmap,
                            show_scalar_bar=False)

            pl.screenshot(f"plots/{modelname}/timeview3D/timeview3D_{type}_{t}.png", 
                          transparent_background=True)
        pl.screenshot(f"plots/{modelname}/timeview3D/timeview3D_{type}.png", 
                      transparent_background=True)


if __name__ == "__main__":
    typer.run(make_3D_time_view)

