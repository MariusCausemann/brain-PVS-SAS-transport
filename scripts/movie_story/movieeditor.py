import os
import typer

from moviepy.editor import (TextClip, ColorClip, CompositeVideoClip,
                            VideoFileClip, concatenate_videoclips) 

# Set background size and color
width = 1920 # pixels
height = 1080 # pixels
black = (0, 0, 0)

def section_clip(title, duration):

    # Create background
    background = ColorClip(size=(width, height), color=black, duration=duration)

    # Create overlay text
    clip = TextClip(title, fontsize=70, color="white", font="Arial-Bold",
                    method="label")
    clip = clip.set_position('center').set_duration(duration)
    clip = clip.fadein(1.0).fadeout(1.0) # 1 sec fade-in and 1 sec fade-out

    # Combine background and text
    composite = CompositeVideoClip([background, clip])
    
    return composite

def generate_full_movie():#, videofile: str):
    
    # A list of clips to be concatenated at the end
    clips = []

    # Main cover: 
    title = "The brain's waterscape: in-silico molecular transport"
    clip = section_clip(title, 5.0)
    clips.append(clip)
    
    # Section #1
    title = "Our subject: healthy male (26 y)"
    clip = section_clip(title, 3.0)
    clips.append(clip)

    # Section #2
    title = "From raw to segmented brain images"
    clip = section_clip(title, 3.0)
    clips.append(clip)
    
    # Section #3
    title = "From segmented images to tetrahedral meshes"
    clip = section_clip(title, 3.0)
    clips.append(clip)

    # Section #4
    title = "From meshes to finite element simulation"
    clip = section_clip(title, 3.0)
    clips.append(clip)

    if False:
        # Add existing video
        # This puts the video on a black background of the correct size
        video_clip = VideoFileClip(videofile)
        video_clip = video_clip.subclip(0, 3) # For more rapid testing
        video_clip.set_position("center")
        video_clip = CompositeVideoClip([video_clip], size=(width, height))
        video_clip = video_clip.fadein(0.5)
        clips.append(video_clip)
    
    # Concatenate everything
    video = concatenate_videoclips(clips)
    
    # Write final file
    os.makedirs("movies", exist_ok=True)
    filename = "movies/insilico_transport.mp4"
    video.write_videofile(filename, fps=24)

    return filename

if __name__ == "__main__":

    typer.run(generate_full_movie)
    
