import pyvista as pv
import numpy as np
import os
import time
import typer


def plot_model(modelname: str, t:int):
    plotdir = f"plots/{modelname}"
    os.makedirs(plotdir, exist_ok=True)

    def getresult(modelname, domain, t):
        filename = f"results/{modelname}/{modelname}_{domain}.pvd"
        reader = pv.get_reader(filename)
        reader.set_active_time_value(t)
        ds = reader.active_datasets[0]
        return pv.read(f"results/{modelname}/{ds.path}")
    
    def time_str(sec):
        m, s = divmod(sec, 60)
        h, m = divmod(m, 60)
        return f"{h:02d}:{m:02d} h"

    sas = getresult(modelname, "sas", t)
    art = getresult(modelname, "arteries", t)
    ven = getresult(modelname, "venes", t)


    clipped = sas.clip(normal="y")


    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(clipped, cmap="matter", clim=(0,1),
                scalar_bar_args=dict(title="conc", vertical=True,
                                    height=0.8,position_y=0.1))
    pl.add_mesh(art, cmap="matter", clim=(0,1), show_scalar_bar=False,
                render_lines_as_tubes=True,line_width=5)
    pl.add_mesh(ven, cmap="matter", clim=(0,1), show_scalar_bar=False,
                render_lines_as_tubes=True,line_width=5)
    pl.camera_position = 'zx'
    pl.camera.roll += 90
    pl.camera.zoom(1.6)

    pl.add_title(f"time: {time_str(t)}", font_size=12)
    pl.screenshot(f"{plotdir}/{modelname}_{t}.png")


if __name__ == "__main__":
    typer.run(plot_model)
