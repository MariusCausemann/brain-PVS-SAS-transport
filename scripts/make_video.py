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

def make_movie(modelname:str):
    config = read_config(f"configfiles/{modelname}.yml")
    dt, T = config["dt"], config["T"]
    #cmap = sns.blend_palette(["lightsteelblue", "orange", "darkviolet"],as_cmap=True)
    cmap = from_k3d(pcm.Yellow___Gray___Blue).reversed()
    clim = (0,2)
    times = np.arange(0, T + dt, dt*config["output_frequency"])
    scalarbar_args = dict(color="white", title="concentration (mmol/l)", 
                          vertical=False, height=0.08, width=0.6,
                          position_x=0.2, position_y=0.05, 
                          title_font_size=36, use_opacity=False,
                          label_font_size=28, fmt="%.1f")
    meshfile = Path(config["mesh"])
    skull = pv.read(meshfile.parent / "surfaces/skull.ply")
    ventricles = pv.read(meshfile.parent / "surfaces/ventricles.ply")
    art = get_result(modelname, "artery", times).ctp()
    sas = get_result(modelname, "sas", times)
    arteries = get_tubes(art)
    arteries["c"] = arteries["c_0"]
    for t in times:
        sas[f"c_{t}"] = np.where(sas[f"c_{t}"] > clim[1], 
                                 clim[1], sas[f"c_{t}"])
    img = pv.create_grid(sas, dimensions=(300, 300, 300))
    img = img.sample(sas, progress_bar=True)

    pl = pv.Plotter(off_screen=True, window_size=(1200, 1200))
    pl.set_background("black")
    pl.open_movie(f"plots/{modelname}/{modelname}.mp4", framerate=24,
                  quality=8)
    pl.add_mesh(skull,style='points',point_size=0.1, opacity=0.6, show_scalar_bar=False,
                color="silver", specular=1, render_points_as_spheres=True)
    pl.add_mesh(ventricles, style='points',point_size=0.8,
                show_scalar_bar=False,specular=1,
                render_points_as_spheres=True,
                color="silver")
    pl.add_mesh(arteries, scalars="c", cmap=cmap, 
                clim=clim, show_scalar_bar=False,)
    pl.camera_position = 'xz'
    pl.camera.zoom = 1.4
    pl.show(auto_close=False) 
    pl.write_frame() 

    # generate frames
    for t in tqdm(list(times[1:]) + [times[-1]]*24*10):
        pl.add_text(f"{int(t/3600)} h", name='time-label', 
                    color="lightgray", font_size=24)
        arteries["c"] = arteries[f"c_{t}"]
        pl.add_volume(img, scalars=f"c_{t}", specular=0.8, diffuse=1, 
                      ambient=0.5,
                      cmap=cmap, opacity=[0, 0.15], clim=clim, shade=False,
                      show_scalar_bar=True, name="vol",
                      opacity_unit_distance=0.0002,
                      scalar_bar_args=scalarbar_args, mapper="gpu")
        pl.camera.azimuth += 360 / len(times) * T/(24*60*60)
        pl.write_frame()  

    pl.close()

if __name__ == "__main__":
    typer.run(make_movie)

