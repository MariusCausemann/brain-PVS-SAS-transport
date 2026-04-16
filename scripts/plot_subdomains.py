import pyvista as pv
import os
from subdomain_ids import PARID, CSFID, LVID,V34ID
import seaborn as sns
import numpy as np

from matplotlib import font_manager

font_path = font_manager.findfont("Nimbus Sans", fallback_to_default=True)

def add_labels_to_scalar_bar(plotter, labels, font_path, font_size=40, y_offset=0.04):
    """
    Automatically centers labels under the plotter's existing scalar bar.
    """
    # 1. Grab the actual VTK actor from the plotter
    sbar = plotter.scalar_bar
    
    # 2. Extract coordinates (relative 0.0 to 1.0)
    # GetPosition() returns [x, y] of bottom-left
    # GetPosition2() returns [width, height]
    pos = sbar.GetPosition()
    size = sbar.GetPosition2()
    
    bar_x = pos[0]
    bar_y = pos[1]
    bar_width = size[0]
    
    n_segments = len(labels)
    segment_width = bar_width / n_segments
    
    for i, label_text in enumerate(labels):
        # Math: Start + (width * percentage across)
        text_x = bar_x + (segment_width * (i + 0.5))
        text_y = bar_y - y_offset
        
        actor = plotter.add_text(label_text, position=(text_x, text_y), 
                                 font_size=font_size,viewport=True)
        
        # Apply the styling
        prop = actor.GetTextProperty()
        prop.SetFontFile(font_path)
        prop.SetJustificationToCentered()
        prop.SetVerticalJustificationToCentered()

def get_camera(grid):
    pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
    pl.add_mesh(grid)
    pl.camera_position = 'yz'
    pl.camera.roll += 0
    pl.camera.azimuth += 30
    #pl.camera.elevation += 10
    pl.camera.zoom(1.4)
    camera =  pl.camera.copy()
    pl.close()
    return camera


if __name__=="__main__":
    grid = pv.read_meshio("mesh/standard/standard.xdmf")
    dirname = "plots/meshplots/"
    os.makedirs(dirname, exist_ok=True)
    camera = get_camera(grid)
    PAR = grid.extract_cells(grid["label"]==PARID)
    CSF = grid.extract_cells(grid["label"]==CSFID).clip()
    V34 = grid.extract_cells(grid["label"]==V34ID)
    LV = grid.extract_cells(grid["label"]==LVID)

    CSF["l"] = np.full(CSF["label"].shape, 0)
    LV["l"] = np.full(LV["label"].shape, 1)
    V34["l"] = np.full(V34["label"].shape, 2)

    comb = pv.merge([CSF, LV, V34])

    bar_args=dict(title=None, vertical=False, height=0.06, width=0.6, position_x=0.2,
                            position_y=0.07, title_font_size=52,
                            bold=False, n_labels=0,
                            label_font_size=76)

    pl = pv.Plotter(off_screen=True, window_size=(1600, 1600))
    #pl.enable_depth_peeling() 
    colors = sns.color_palette("Blues", n_colors=3)

    pl.add_mesh(comb, scalars="l", show_scalar_bar=True, 
                             scalar_bar_args=bar_args,
                             cmap="Blues",
                categories=True,
                #annotations={0:"CSF", 1:"LV", 2:"V3&4"}
                )
    
    pl.add_mesh(grid.extract_surface(), opacity=0.4, color=colors[0],
                show_scalar_bar=False, scalars=None, reset_camera=False)
    
    labels = ["CSF", "LV", "V3&4"]
    add_labels_to_scalar_bar(pl, labels, font_path, font_size=46)

    pl.camera = camera
    pl.screenshot(f"{dirname}/subdomains.png",
    transparent_background=True)
