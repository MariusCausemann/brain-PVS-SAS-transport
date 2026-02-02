# Author: Marie E. Rognes (meg@simula.no)
# License: MIT
#
# Make sure to have dependencies installed/available:
# conda activate pvs_transport_env 
# 
#
# How to run: 
# python3 scripts/movie_story/make_story_video.py .
#

#import pyvista as pv
import typer
from pathlib import Path

from moviepy.editor import TextClip, ColorClip, CompositeVideoClip

def make_movie(framedir: str):
    # This overview function simply collects all frames (images),
    # connects them with text clips, and uses moviepy to 
    
    # Add all clips to this and composite video at the end
    clips = []

    # Set background size and color
    width = 1920 # pixels
    height = 1080 # pixels
    black = (0, 0, 0)
    
    # Main cover
    duration = 5.0
    background = ColorClip(size=(width, height), color=black, duration=duration)

    # Main cover 
    title = "The brain's waterscape: in-silico molecular transport"
    clip = TextClip(title, fontsize=70, color="white", font="Arial-Bold",
                    method="label")
    clip = clip.set_position('center').set_duration(5)
    clip = clip.fadein(1.0) # This adds a 1-second fade-in from black 
    clips.extend([background, clip])
    
    # Composite everything
    video = CompositeVideoClip(clips)

    # Write final file
    filename = "insilico_transport.mp4"
    video.write_videofile(filename, fps=24)

    return filename

if __name__ == "__main__":

    typer.run(make_movie)

    
# def add_heading(pl, heading, frames):
#     # Show the heading for the given number of frames 
#     pl.add_text(heading, position=(200, 600), name="heading", font_size=24)
#     pl.show(auto_close=False)
#     for i in range(frames):
#         pl.write_frame()  

# def add_T1w_images(pl, datadir):

#     # Open and read the data first 
#     filename = Path(datadir, "data", "T1_synthseg.nii.gz")
    
#     print("Reading image from ", filename)
#     reader = pv.get_reader(filename) 
#     img = reader.read()
#     print(img)

    # # Set global plotting options for Pyvista
    # pv.set_plot_theme("dark")
    # pv.global_theme.font.family = "arial"
    # pv.global_theme.font.color = "white"

