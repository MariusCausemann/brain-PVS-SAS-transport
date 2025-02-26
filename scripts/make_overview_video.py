# Make sure that dependencies are in place first:
#
# $ conda activate pvs_transport_env
#
from plotting_utils import read_config
import pyvista as pv
import typer
from pathlib import Path
import numpy as np
# from tqdm import tqdm
# from plot_csf_flow import from_k3d
import k3d.colormaps.paraview_color_maps as pcm
# from extract_vessels import get_tubes
# import seaborn as sns

def show_heading(pl, text, frames, clear=True):

    # Show the heading for 'frames' frames
    pl.add_text(text, position=(200, 600), name="heading", font_size=24)
    pl.show(auto_close=False)
    for i in range(frames):
        pl.write_frame()  

    # Clear the text
    pl.add_text("", position=(200, 600), name="heading", font_size=24)
    pl.show(auto_close=False)
    for i in range(2):
        pl.write_frame()  
    
def make_overview_movie(datadir: str):

    if False:
        config = read_config(f"configfiles/modelA.yml")
        dt, T = config["dt"], config["T"]
        times = np.arange(0, T + dt, dt*config["output_frequency"])
        print("T:", T, "\ntimes (s):", times)
    
    pv.set_plot_theme("dark")
    pv.global_theme.font.family = "arial"
    #pv.global_theme.font.color = "white"

    # Show the original T1w MRI images (almost)
    filename = Path(datadir, "data", "T1_synthseg_robust.nii.gz")
    print("Reading image from ", filename)
    reader = pv.get_reader(filename) 
    img = reader.read()

    # Open the movie
    pl = pv.Plotter(off_screen=True, window_size=(1920, 1200))
    fps = 24
    pl.open_movie("overview.mp4", framerate=fps, quality=8)

    # Subject details as intro
    show_heading(pl, "Healthy male (26 y)", 2*fps)

    text = "Magnetic resonance images (3T, T1w, 0.8mm x 0.8mm x 0.8mm)"
    show_heading(pl, text, 2*fps)
    
    slices = img.slice_orthogonal()
    pl.add_mesh(slices, cmap="gist_gray", show_scalar_bar=False)
    #pl.add_mesh(slices, cmap="tab20c")
    #pl.add_volume(img, cmap="tab20c")#, show_scalar_bar=True)
    pl.camera_position = 'iso'
    pl.camera.zoom(1.6)
    pl.show(auto_close=False)
    for i in range(fps):
        pl.write_frame()  

    # Rotate image to see details
    rotations = 2*fps
    for i in range(rotations):
        pl.camera.azimuth = i*360 / rotations #* T/(24*60*60)
        print("Setting camera azimuth to", pl.camera.azimuth)
        #pl.write_frame()  
        pl.show(auto_close=False)
        #pl.screenshot("foo%d.png" % i)
        pl.write_frame()  
        pl.write_frame()  

    pl.close()
    exit()
    
    #pl.open_movie("foo.mp4", quality=10)
    #pl.add_volume(img, cmap="bone", show_scalar_bar=False)


    for i in range(22):

        pl = pv.Plotter(off_screen=True, window_size=(1920, 1200))
        #pl.set_background('black')
        
        #pl.open_movie("foo.mp4", quality=10)
        #pl.add_volume(img, cmap="bone", show_scalar_bar=False)
        
        print("i = ", i)
        slices = img.slice_orthogonal(y=9*i)
        pl.add_text("%s" % i, name='time-label', 
                    color="lightgray", font_size=24)
        pl.add_mesh(slices, cmap="bone", show_scalar_bar=False)

        pl.camera.zoom(1)
        pl.show(auto_close=False)
        #pl.write_frame() 
        pl.screenshot("foo_%d.png" % i)

    pl.close()
    print("Video done")
    
if __name__ == "__main__":
    typer.run(make_overview_movie)
