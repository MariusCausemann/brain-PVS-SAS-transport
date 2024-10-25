import pyvista as pv
import os
from generate_synthseg_mesh import PARID, CSFID, LVID,V34ID
import seaborn as sns
import numpy as np

grid = pv.read_meshio("mesh/standard/standard.xdmf")
dirname = "plots/meshplots/"
os.makedirs(dirname, exist_ok=True)

PAR = grid.extract_cells(grid["label"]==PARID)
CSF = grid.extract_cells(grid["label"]==CSFID)
V34 = grid.extract_cells(grid["label"]==V34ID)
LV = grid.extract_cells(grid["label"]==LVID)

CSF["l"] = np.full(CSF["label"].shape, " SAS", dtype='<U10')
V34["l"] = np.full(V34["label"].shape, "V3 & V4", dtype='<U10')
LV["l"] = np.full(LV["label"].shape, "LV", dtype='<U10')

comb = pv.merge([CSF.clip(), V34, LV])

bar_args=dict(title="", vertical=False, height=0.08, width=0.6, position_x=0.2,
                        position_y=-0.05, title_font_size=52,
                        bold=False, font_family="times",
                        label_font_size=44)

pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
pl.add_mesh(CSF.clip())
pl.camera_position = 'yz'
pl.camera.roll += 0
pl.camera.azimuth += 30
#pl.camera.elevation += 10
pl.camera.zoom(1.2)
camera = pl.camera.copy()

pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
colors = sns.color_palette("Blues", n_colors=3)
pl.add_mesh(comb, scalars="l", show_scalar_bar=True, cmap="Blues",
            scalar_bar_args=bar_args)
pl.scalar_bar.SetTextPosition(0)
pl.camera = camera
pl.screenshot(f"{dirname}/subdomains.png",
transparent_background=True)
