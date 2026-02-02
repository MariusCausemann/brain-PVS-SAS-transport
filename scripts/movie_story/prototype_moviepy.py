# conda activate pvs_transport_env

from moviepy.editor import TextClip, ColorClip, CompositeVideoClip

clips = []

# Set background size and color
width = 1920 # pixels
height = 1080 # pixels
black = (0, 0, 0)

# Define the background
duration = 5.0
background = ColorClip(size=(width, height), color=black, duration=duration) 

# Create a text overlay 
title = TextClip("Title", fontsize=70, color="white", font="Arial-Bold",
                 method="label")
title = title.set_position('center').set_duration(duration)
title = title.fadein(1.0) # This adds a 1-second fade-in from black 
clips.extend([background, title])

# Add all clips together into a composite video
video = CompositeVideoClip(clips)

# Write final file
video.write_videofile("moviepy_prototype.mp4", fps=24)
