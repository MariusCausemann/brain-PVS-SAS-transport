import pyvista as pv
import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def mycolormap(twilight=False):
    if twilight:
        original_cmap = plt.get_cmap("twilight_shifted")
        colors = original_cmap(np.linspace(0, 1, 256))
        
        # Force the first color (index 0) to be pure black
        colors[0] = [0, 0, 0, 1]

        # Create a new matplotlib colormap object
        colormap = ListedColormap(colors)
        name = "twilight"
    else:
        colormap = "gist_gray"
        name = colormap
        
    return (colormap, name)


#filename = "data/T1.nii.gz"
filename = "data/T1_synthseg.nii.gz"
#(colormap, name) = mycolormap("gist_gray")
(colormap, name) = mycolormap("twilight")

width = 1920
height = 1080

# Setup folders
output_dir = f"frames/mri_T1_synthseg_slices_{name}"
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

for k in range(1, z_max, 1):
    
    # Extract single slice
    print("k = ", k)
    slice = img.slice(normal=(0, 0, 1), origin=(0, 0, k))

    plotter.clear()
    plotter.add_mesh(outline, opacity=0.0) # Keep camera locked
    
    # Render slice
    plotter.add_mesh(slice, cmap=colormap,
                     show_scalar_bar=False, lighting=False)
    
    # Save to output_dir/
    plotter.screenshot(f"{output_dir}/frame_{k:04d}.png")
    
print(f"Done! Check {output_dir}.")
