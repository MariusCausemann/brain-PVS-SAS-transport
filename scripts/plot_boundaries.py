import pyvista as pv
import os
from mark_and_refine_mesh import (PIA_ID, LV_INTERF_ID, SPINAL_OUTLET_ID,
LOWER_SKULL_ID, UPPER_SKULL_ID, SPINAL_CORD_ID)
import seaborn as sns
import numpy as np
from cmap import Colormap
from plot_subdomains import get_camera

grid = pv.read_meshio("mesh/standard/standard_facets.xdmf")
dirname = "plots/meshplots/"
os.makedirs(dirname, exist_ok=True)

bds = {"AM-U":UPPER_SKULL_ID,"AM-L": LOWER_SKULL_ID,
       "SSAS":SPINAL_OUTLET_ID,
        #"SC": SPINAL_CORD_ID, 
        " LV":LV_INTERF_ID,  " Pia":PIA_ID,}
clips = [" Pia", "AM-L", "AM-U"]

bar_args=dict(title="", vertical=True, height=0.8, width=0.06, position_x=-0.035,
                        position_y=0.1, title_font_size=52,
                        bold=False, font_family="times",
                        label_font_size=72)
print(grid.array_names)
bounds = list(grid.bounds)
bounds[0] = (bounds[0] + bounds[1]) * 0.5
bounds[2] = (bounds[2] + bounds[3]) * 0.4
print(bounds)
colors = sns.color_palette("twilight_shifted", n_colors=20, as_cmap=True)
#colors = sns.diverging_palette(145, 300, s=60, as_cmap=True)

meshes= []
for i, (n, bid) in enumerate(bds.items()):
    m = grid.extract_cells(grid["f"]==bid)
    if n in clips:
        m = m.clip() #m.clip_box(bounds)
    m["l"] = np.full(m["f"].shape, n, dtype='<U10')
    meshes.append(m)


camera = get_camera(grid)

pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
pl.add_mesh(pv.merge(meshes), scalars="l", show_scalar_bar=True, 
            scalar_bar_args=bar_args, 
            cmap=Colormap('cmocean:curl_pink_r'),
            #cmap=Colormap('vispy:pugr')
            )
pl.scalar_bar.SetTextPosition(0)
pl.camera = camera
pl.screenshot(f"{dirname}/boundaries.png",
transparent_background=True)
print(pl.camera.position)