import pyvista as pv
import os
from mark_and_refine_mesh import (PIA_ID, LV_INTERF_ID, SPINAL_OUTLET_ID,
LOWER_SKULL_ID, UPPER_SKULL_ID, SPINAL_CORD_ID)
import seaborn as sns
import numpy as np
from cmap import Colormap
from plot_subdomains import get_camera
from matplotlib import font_manager
from plot_subdomains import add_labels_to_scalar_bar

font_path = font_manager.findfont("Nimbus Sans", fallback_to_default=False)

grid = pv.read_meshio("mesh/standard/standard_facets.xdmf")
dirname = "plots/meshplots/"
os.makedirs(dirname, exist_ok=True)

bds = {"AM-U":UPPER_SKULL_ID,"AM-L": LOWER_SKULL_ID,
       "SSAS":SPINAL_OUTLET_ID,
        #"SC": SPINAL_CORD_ID, 
        "LV":LV_INTERF_ID,  "Pia":PIA_ID,}
clips = ["Pia", "AM-L", "AM-U"]

bar_args=dict(title="", vertical=False, height=0.06, width=0.85, position_x=0.05,
                        position_y=0.07, title_font_size=52,
                        bold=False, n_labels=0,
                        label_font_size=76)
print(grid.array_names)
bounds = list(grid.bounds)
bounds[0] = (bounds[0] + bounds[1]) * 0.5
bounds[2] = (bounds[2] + bounds[3]) * 0.4
print(bounds)
colors = sns.color_palette("twilight_shifted", n_colors=20, as_cmap=True)
#colors = sns.diverging_palette(145, 300, s=60, as_cmap=True)

cat_map = {"LV": 0, "Pia": 1, "AM-L": 2, "AM-U": 3, "SSAS": 4}
meshes= []
for i, (n, bid) in enumerate(bds.items()):
    m = grid.extract_cells(grid["f"]==bid)
    if n in clips:
        m = m.clip() #m.clip_box(bounds)
    m["l"] = np.full(m["f"].shape, cat_map[n])
    meshes.append(m)
#from IPython import embed; embed()
camera = get_camera(grid)

pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
pl.add_mesh(pv.merge(meshes), scalars="l", show_scalar_bar=True, 
            scalar_bar_args=bar_args, 
            cmap=Colormap('cmocean:curl_pink_r').to_matplotlib(),
            categories=True,
            )
add_labels_to_scalar_bar(pl, cat_map.keys(), font_path, font_size=46)
pl.camera = camera
pl.screenshot(f"{dirname}/boundaries.png",
transparent_background=True)
print(pl.camera.position)