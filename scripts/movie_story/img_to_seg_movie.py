from moviepy.editor import ImageSequenceClip, CompositeVideoClip, TextClip
from moviepy.editor import concatenate_videoclips

# 1. Load the raw slices
print("Loading frames/mri_T1_slices_gist_gray")
clip_raw = ImageSequenceClip("frames/mri_T1_slices_gist_gray", fps=24)

# Just take part of it, fade in and fade out
duration = clip_raw.duration
clip_raw = clip_raw.subclip(2.0, clip_raw.duration-0.5)
clip_raw = clip_raw.fadein(0.5).fadeout(0.5)

# Add overlay text
text_clip = TextClip("T1-weighted MRI\n (1 x 1 x 1 mm^3)",
                     fontsize=60, color="white", font="Helvetica",
                     size=(600, None),
                     method="caption")
text_clip = text_clip.set_position((50, 800))
text_clip = text_clip.set_duration(clip_raw.duration)
text_clip = text_clip.fadein(1.0).fadeout(1.0) # 2 sec fade-in, 1 sec fade-out

# Combine background and text
part1 = CompositeVideoClip([clip_raw, text_clip])

# 2. Load the segmented images
print("Loading frames/mri_T1_synthseg_slices_twilight")
clip = ImageSequenceClip("frames/mri_T1_synthseg_slices_twilight", fps=24)
clip = clip.fadein(1.0).fadeout(1.0)

# Add overlay text
text_clip = TextClip("Segmentation into brain regions via SynthSeg",
                     fontsize=60, color="white", font="Helvetica",
                     size=(600, None),
                     method="caption")
text_clip = text_clip.set_position((1200, "center"))
text_clip = text_clip.set_duration(clip.duration)
text_clip = text_clip.fadein(1.0).fadeout(0.5) # 2 sec fade-in, 1 sec fade-out

# Combine background and text
part2 = CompositeVideoClip([clip, text_clip])

#final_video = concatenate_videoclips([part1, part2])#, part2, part3])
final_video = concatenate_videoclips([part1, part2])

final_video.write_videofile("01_images_to_seg.mp4", fps=24)
