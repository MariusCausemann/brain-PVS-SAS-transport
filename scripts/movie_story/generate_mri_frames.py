import pyvista as pv
import numpy as np
import os

filename = "data/T1_synthseg.nii.gz"
width = 1920
height = 1080

# Setup folders
colormap = "twilight_shifted"
output_dir = f"frames/mri_frames_{colormap}"
os.makedirs(output_dir, exist_ok=True)

# 1. Load data
img = pv.read(filename)

# 2. Setup plotter
plotter = pv.Plotter(off_screen=True, window_size=[1920, 1080])
plotter.set_background("black")

# Add an invisible outline box to lock the camera zoom/position
outline = img.outline()
plotter.add_mesh(outline, opacity=0.0) 
plotter.camera_position = 'iso'
plotter.camera.zoom(1.5)

# 3. Animation loop
dims = img.dimensions
spacing = img.spacing
origin = img.origin
z_max = dims[2]

for k in range(0, z_max, 10):

    print("k", k)
    # We slice the raw MRI volume
    slice = img.slice_orthogonal(z=k)
    
    plotter.clear()
    plotter.add_mesh(outline, opacity=0.0) # Keep camera locked
    
    # Render raw slice (Grayscale)
    plotter.add_mesh(slice, cmap=colormap,
                     show_scalar_bar=False, lighting=False)
    
    # Save to output_dir/
    plotter.screenshot(f"{output_dir}/frame_{k:04d}.png")
    
print(f"Done! Check {output_dir}.")
