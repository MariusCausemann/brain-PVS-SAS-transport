from moviepy.editor import ImageSequenceClip, CompositeVideoClip, TextClip
from moviepy.editor import concatenate_videoclips

# 1. Load the raw slices
clip_raw = ImageSequenceClip("frames/mri_slices_T1raw_gist_gray", fps=24)

# Just take part of it, fade in and fade out
duration = clip_raw.duration
clip_raw = clip_raw.subclip(2.0, clip_raw.duration-0.5).fadein(1.0).fadeout(1.0)

# Add overlay text
text_clip = TextClip("T1-weighted MRI\n (1 x 1 x 1 mm^3)",
                     fontsize=60, color="white", font="Helvetica",
                     size=(600, None),
                     method="caption")
text_clip = text_clip.set_position((50, 800))
text_clip = text_clip.set_duration(clip_raw.duration)
text_clip = text_clip.fadein(2.0).fadeout(1.0) # 1 sec fade-in and 1 sec fade-out

# Combine background and text
part1 = CompositeVideoClip([clip_raw, text_clip])

# 2. Define the Transition Timing
# Let's say the full scan is 6 seconds.
# We want to switch from Gray to Color at 3 seconds.
#fade_start = 4.0
#fade_duration = 1.0 

# # 3. Composite
# # We overlay the Color clip on top of the Raw clip.
# # We tell the Color clip to be invisible (opacity 0) until 'fade_start',
# # then fade in to opacity 1.0 over 'fade_duration'.

# clip_seg_fading = clip_seg.set_start(0).crossfadein(fade_duration)
# # Note: Complex crossfades in MoviePy can be tricky. 
# # A simpler, robust way is to slice them:

# # --- ROBUST METHOD ---
# # Part 1: Pure Raw (0 to 2.5s)

# # Part 2: The Crossfade (2.5s to 3.5s)
# # We take the 1-second chunk from both and blend them
# part2_raw = clip_raw.subclip(fade_start, fade_start + fade_duration)
# part2_seg = clip_seg.subclip(fade_start, fade_start + fade_duration)
# part2 = CompositeVideoClip([part2_raw, part2_seg.crossfadein(fade_duration)])

# # Part 3: Pure Seg (3.5s to End)
# part3 = clip_seg.subclip(fade_start + fade_duration)

# Stitch them
final_video = concatenate_videoclips([part1])#, part2, part3])

final_video.write_videofile("01_images_to_seg_raw.mp4", fps=24)
