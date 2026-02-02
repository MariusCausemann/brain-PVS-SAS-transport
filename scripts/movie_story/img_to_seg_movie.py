from moviepy.editor import ImageSequenceClip, CompositeVideoClip

# 1. Load the two sequences
# fps=30 makes the scan smooth
clip_raw = ImageSequenceClip("frames/mri_frames_gist_gray", fps=2)
clip_seg = ImageSequenceClip("frames/mri_frames_twilight_shifted", fps=2)

# 2. Define the Transition Timing
# Let's say the full scan is 6 seconds.
# We want to switch from Gray to Color at 3 seconds.
fade_start = 4.0
fade_duration = 1.0 

# 3. Composite
# We overlay the Color clip on top of the Raw clip.
# We tell the Color clip to be invisible (opacity 0) until 'fade_start',
# then fade in to opacity 1.0 over 'fade_duration'.

clip_seg_fading = clip_seg.set_start(0).crossfadein(fade_duration)
# Note: Complex crossfades in MoviePy can be tricky. 
# A simpler, robust way is to slice them:

# --- ROBUST METHOD ---
# Part 1: Pure Raw (0 to 2.5s)
part1 = clip_raw.subclip(0, fade_start)

# Part 2: The Crossfade (2.5s to 3.5s)
# We take the 1-second chunk from both and blend them
part2_raw = clip_raw.subclip(fade_start, fade_start + fade_duration)
part2_seg = clip_seg.subclip(fade_start, fade_start + fade_duration)
part2 = CompositeVideoClip([part2_raw, part2_seg.crossfadein(fade_duration)])

# Part 3: Pure Seg (3.5s to End)
part3 = clip_seg.subclip(fade_start + fade_duration)

# Stitch them
from moviepy.editor import concatenate_videoclips
final_video = concatenate_videoclips([part1, part2, part3])

final_video.write_videofile("01_scan_transition.mp4", fps=30)
