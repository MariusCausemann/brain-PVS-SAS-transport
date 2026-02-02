# Author: Marie E. Rognes (meg@simula.no)
# License: MIT
#
# Make sure to have dependencies installed/available:
# conda activate pvs_transport_env 
#
# How to run: 
# python3 scripts/make_story_video.py .
#

import pyvista as pv
import typer
from pathlib import Path

def add_heading(pl, heading, frames):
    # Show the heading for the given number of frames 
    pl.add_text(heading, position=(200, 600), name="heading", font_size=24)
    pl.show(auto_close=False)
    for i in range(frames):
        pl.write_frame()  

def add_T1w_images(pl, datadir):

    # Open and read the data first 
    filename = Path(datadir, "data", "T1_synthseg.nii.gz")
    
    print("Reading image from ", filename)
    reader = pv.get_reader(filename) 
    img = reader.read()
    print(img)

def make_movie(imagedir: str):

    # Initialize the plotter and movie
    pl = pv.Plotter(off_screen=True, window_size=(1920, 1200))
    fps = 24 #? 
    qlty = 10 #?
    filename = "personalized_images.mp4"
    pl.open_movie(filename, framerate=fps, quality=qlty)

    # Subject details as intro
    add_heading(pl, "Subject: healthy male (26 y)", 2*fps)

    add_T1w_images(pl, imagedir)
    
    # Wrap up
    pl.close()
    print(filename, "done")

    return filename
    
if __name__ == "__main__":

    # Set global plotting options for Pyvista
    pv.set_plot_theme("dark")
    pv.global_theme.font.family = "arial"
    pv.global_theme.font.color = "white"

    typer.run(make_movie)
