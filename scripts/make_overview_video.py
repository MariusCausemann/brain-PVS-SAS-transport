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
#import seaborn as sns

def show_heading(pl, heading, frames, subtitle=None):

    # Show the heading for 'frames' frames
    pl.add_text(heading, position=(200, 600), name="heading", font_size=24)
    if subtitle:
        pl.add_text(subtitle, position=(200, 400),name="subtitle", font_size=24)
    pl.show(auto_close=False)
    for i in range(frames):
        pl.write_frame()  

    # Clear the text
    pl.add_text("", position=(200, 600), name="heading", font_size=24)
    if subtitle:
        pl.add_text("", position=(200, 400), name="subtitle", font_size=24)
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
    filename = Path(datadir, "data", "T1_synthseg.nii.gz")
    print("Reading image from ", filename)
    reader = pv.get_reader(filename) 
    img = reader.read()
    print(img)
    
    # Open the movie
    pl = pv.Plotter(off_screen=True, window_size=(1920, 1200))
    fps = 24
    pl.open_movie("overview.mp4", framerate=fps, quality=10)

    # Subject details as intro
    show_heading(pl, "Subject: healthy male (26 y)", 2*fps)

    if False:
        heading = "Magnetic resonance images (T1-w)"
        #subtitle = "(3T, 0.8mm x 0.8mm x 0.8mm)"
        show_heading(pl, heading, 2*fps)
    
        slices = img.slice_orthogonal()
        pl.add_mesh(slices, cmap="gist_gray", show_scalar_bar=False)
        pl.camera_position = 'iso'
        pl.camera.zoom(1.6)
        pl.show(auto_close=False)
        for i in range(fps):
            pl.write_frame()  

        # Rotate image to see details
        rotations = 8*fps
        for i in range(rotations):
            pl.camera.azimuth = i*360 / rotations #* T/(24*60*60)
            print("Setting camera azimuth to", pl.camera.azimuth)
            pl.show(auto_close=False)
            pl.write_frame()  

        for i in range(fps):
            pl.write_frame()  

        # Pause a bit
        pl.clear()
        for j in range(2):
            pl.write_frame()  

    if False:
        heading = "Segment brain regions and CSF spaces"
        show_heading(pl, heading, 2*fps)

        # Show segmentation more clearly
        for y in range(20, 180, 1):
            print("y = ", y)
            slices = img.slice_orthogonal(y=y)
            pl.clear()
            pl.add_mesh(slices, cmap="twilight_shifted", lighting=True,
                        show_scalar_bar=False)
            pl.camera_position = 'iso'
            pl.camera.zoom(1.6)
            pl.show(auto_close=False)
            pl.write_frame()  


    # Show mesh 
    pl.clear()
    for j in range(2):
        pl.write_frame()  

    # Plot surfaces
    LV = pv.read(Path(datadir, "mesh/standard/surfaces/LV.ply"))
    V34 = pv.read(Path(datadir,  "mesh/standard/surfaces/V34.ply"))
    parenchyma = pv.read(Path(datadir, "mesh/standard/surfaces/parenchyma.ply"))
    skull = pv.read(Path(datadir, "mesh/standard/surfaces/skull.ply"))

    pl.add_mesh(skull, show_scalar_bar=False, opacity=0.01) # Just for scale
    pl.add_mesh(V34, show_scalar_bar=False, opacity=0.6, color="blue")
    
    pl.camera.zoom(2.0)
    pl.show(auto_close=False)

    pl.add_mesh(LV, show_scalar_bar=False, opacity=0.6, color="blue")

    viewup = [0, 0, 1]
    path = pl.generate_orbital_path(factor=2.0, n_points=8*fps,
                                    viewup=viewup, shift=0.2)
    # Do e.g. half path here., see PolyData doc.
    print(path.points)
    print(path.lines)
    exit()
    pl.orbit_on_path(path, write_frames=True, viewup=viewup, progress_bar=True)
    
    pl.close()
    print("Video done")
    
if __name__ == "__main__":
    typer.run(make_overview_movie)
