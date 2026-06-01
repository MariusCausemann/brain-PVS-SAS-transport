from plotting_utils import get_result, read_config
import pyvista as pv
import typer
import numpy as np
from pathlib import Path
from tqdm import tqdm
from plot_csf_flow import from_k3d
import k3d.colormaps.paraview_color_maps as pcm
from extract_vessels import get_tubes
import numpy as np

def time_interpolate(grid, scalar_str, t_vid, sim_dt, T):
    # Calculate weights based on simulation time
    idx_a = int(t_vid // sim_dt)
    idx_b = idx_a + 1
        
    t_a, t_b = idx_a * sim_dt, idx_b * sim_dt 
    if t_b > T: t_b = t_a   
    w = (t_vid - t_a) / (t_b - t_a) if t_b > t_a else 0

    # Interpolate using the loaded subset keys
    return (1 - w) * grid[scalar_str.format(t=t_a)] + w * grid[scalar_str.format(t=t_b)]


def make_movie(modelname:str, target_dt: float = 300.0):
    config = read_config(f"configfiles/{modelname}.yml")
    sim_dt , T= config["dt"]*config["output_frequency"], config["T"]
    #cmap = sns.blend_palette(["lightsteelblue", "orange", "darkviolet"],as_cmap=True)
    video_times = np.arange(0, T + target_dt, target_dt)
    indices_a = (video_times // sim_dt).astype(int)

    # 2. The 'ceiling' is just the next step
    indices_b = indices_a + 1

    # 3. Combine them and keep only unique, valid indices
    max_idx = int(T // sim_dt)
    required_indices = np.unique(np.concatenate([indices_a, indices_b]))
    required_indices = required_indices[required_indices <= max_idx]

    # 4. Convert back to time values for your get_result call
    needed_sim_times = (required_indices * sim_dt).tolist()
    cmap = from_k3d(pcm.Yellow___Gray___Blue).reversed()

    clim = (0,2)
    N = 300
    scalarbar_args = dict(color="white", title="concentration (mmol/l)", 
                          vertical=False, height=0.08, width=0.6,
                          position_x=0.2, position_y=0.05, 
                          title_font_size=36, use_opacity=False,
                          label_font_size=28, fmt="%.1f")
    meshfile = Path(config["mesh"])
    skull = pv.read(meshfile.parent / "surfaces/skull.ply")
    ventricles = pv.read(meshfile.parent / "surfaces/ventricles.ply")
    art = get_result(modelname, "artery", needed_sim_times).ctp()
    sas = get_result(modelname, "sas", needed_sim_times)
    arteries = get_tubes(art)
    arteries["c"] = arteries["c_0"]
    for t in needed_sim_times:
        sas[f"c_{t}"] = np.where(sas[f"c_{t}"] > clim[1], 
                                 clim[1], sas[f"c_{t}"])
    img = pv.create_grid(sas, dimensions=(N, N, N))
    img = img.sample(sas, progress_bar=True)

    pl = pv.Plotter(off_screen=True, window_size=(1200, 1200))
    pl.set_background("black")
    pl.open_movie(f"plots/{modelname}/{modelname}.mp4", framerate=24, quality=8)
    pl.add_mesh(skull,style='points',point_size=0.1, opacity=0.6, show_scalar_bar=False,
                color="silver", specular=1, render_points_as_spheres=True)
    pl.add_mesh(ventricles, style='points',point_size=0.8,
                show_scalar_bar=False,specular=1,
                render_points_as_spheres=True,
                color="silver")
    pl.add_mesh(arteries, scalars="c", cmap=cmap, 
                clim=clim, show_scalar_bar=False,)
    pl.camera_position = 'xz'
    pl.camera.focal_point = np.array(skull.center) - np.array([0,0,0.011])
    pl.camera.zoom(1.4)

    img["vol_scalars"] = img["c_0"]
    volume = pl.add_volume(img, scalars=f"vol_scalars", specular=0.8, diffuse=1, 
                    ambient=0.5,
                    cmap=cmap, opacity=[0, 0.15], clim=clim, shade=False,
                    show_scalar_bar=True, name="vol",
                    opacity_unit_distance=0.0002,
                    scalar_bar_args=scalarbar_args, mapper="gpu")
    volume.mapper.SetAutoAdjustSampleDistances(True)
    
    pl.show(auto_close=False) 
    pl.write_frame() 


    # generate frames
    for t in tqdm(list(video_times[1:]) + [video_times[-1]]*24*10):
        pl.add_text(f"{int(t/3600)} h", name='time-label', 
                    color="lightgray", font_size=24)
        arteries["c"][:] = time_interpolate(arteries, "c_{t}", t, sim_dt, T)
        img["vol_scalars"][:] = time_interpolate(img, "c_{t}", t, sim_dt, T)
        pl.camera.azimuth += 360 / len(video_times) * T/(24*60*60)
        pl.write_frame()  

    pl.close()

if __name__ == "__main__":
    typer.run(make_movie)

