import pyvista as pv
import yaml
import collections
import matplotlib.pyplot as plt

def set_plotting_defaults():
    import seaborn as sns
    plt.rc("axes.spines", top=False, right=False)
    sns.set_context("talk")
    sns.set_palette("blend:#7AB,#EDA")

def get_result(modelname, domain, times):
    filename = f"results/{modelname}/{modelname}_{domain}.pvd"
    reader = pv.get_reader(filename)
    if not isinstance(times, collections.Iterable):
        times = [times]
    reader.set_active_time_value(times[0])
    ds = reader.active_datasets[0]
    data = pv.read(f"results/{modelname}/{ds.path}")
    if len(times) > 1:
        for ar in data.array_names:
            data.rename_array(ar, f"{ar}_{times[0]}")
    for t in times:
        reader.set_active_time_value(t)
        ds = reader.active_datasets[0]
        d =  pv.read(f"results/{modelname}/{ds.path}")
        for ar in d.array_names:
            data[f"{ar}_{t}"] = d[ar]
    return data

def time_str(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return f"{h:02d}:{m:02d}"

def read_config(configfile):
    with open(configfile) as conf_file:
        config = yaml.load(conf_file, Loader=yaml.FullLoader)
    return config


def clip_plot(sas, networks, filename, title, clim, cmap, cbar_title):
    clipped = sas.clip(normal="y")
    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(clipped, cmap=cmap, clim=clim,
                scalar_bar_args=dict(title=cbar_title, vertical=False,
                                    height=0.1, width=0.6, position_x=0.2,
                                    position_y=0.0, title_font_size=36,
                                    label_font_size=32))
    for netw in networks:
        pl.add_mesh(netw, cmap=cmap, clim=clim, show_scalar_bar=False,
                    render_lines_as_tubes=True,line_width=5)

    pl.camera_position = 'zx'
    pl.camera.roll += 90
    pl.camera.zoom(1.6)
    pl.add_title(title, font_size=12)
    pl.screenshot(filename, transparent_background=True)
