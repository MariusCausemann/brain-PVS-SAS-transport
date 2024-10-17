import numpy as np
import yaml
import collections
from typing import List
from IPython import embed

def set_plotting_defaults():
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.rc("axes.spines", top=False, right=False)
    sns.set_context("notebook")
    sns.set_palette("blend:#7AB,#EDA")

def minmax(arr_list, percentile=95):
    if percentile is None:
        return float(np.min(arr_list)), float(np.max(arr_list))
    else:
        percentiles = [np.percentile(arr, [100 - percentile, percentile])
                       for arr in arr_list]
        return float(np.min(percentiles)), float(np.max(percentiles))

def get_result(modelname, domain, times):
    from vtk_adapter import create_vtk_structures
    import pyvista as pv
    if not isinstance(times, collections.Iterable):
        times = [times]
    res, label = get_result_fenics(modelname, domain, times)
    Vh = res[0].function_space()
    topology, cell_types, x = create_vtk_structures(Vh)
    grid = pv.UnstructuredGrid(topology, cell_types, x)
    for r,t in zip(res, times):
        grid[f"{r.name()}_{t}"] = r.vector()[:]
    if domain == "sas":
        grid["label"] = label.array()[:]
    else:
        grid["radius"] = label.vector()[:]
    return grid

def get_result_fenics(modelname, domain, times):
    from fenics import XDMFFile, FunctionSpace, Function
    from solver  import get_mesh, read_vtk_network, as_P0_function
    filename = f"results/{modelname}/{modelname}_{domain}.xdmf"
    config = read_config(f"configfiles/{modelname}.yml")
    dt , T= config["dt"], config["T"]
    alltimes = np.arange(0, T + dt, dt*config["output_frequency"])
    if domain=="sas":
        mesh, vol_subdomains = get_mesh(config["mesh"])
        V = FunctionSpace(mesh, "DG", 1)
    if domain == "artery":
        mesh, radii, _ = read_vtk_network("mesh/networks/arteries_smooth.vtk", rescale_mm2m=False)
        V = FunctionSpace(mesh, "CG", 1)
    if domain == "vein":
        mesh, radii, _ = read_vtk_network("mesh/networks/venes_smooth.vtk", rescale_mm2m=False)
        V = FunctionSpace(mesh, "CG", 1)
    results = []
    with XDMFFile(filename) as f:
        for t in times:
            i = np.where(alltimes==t)[0][0]
            c = Function(V)
            c.rename("c","c")
            f.read_checkpoint(c, "c", i)
            results.append(c)
    if domain in ["artery", "vein"]:
        return results, as_P0_function(radii)
    elif domain=="sas":
        return results, vol_subdomains

def get_range(modelname, domain, times, var, percentile=95):
    grid = get_result(modelname, domain, times)
    return minmax([grid[f"{var}_{t}"] for t in times], percentile=percentile)

def get_diff_range(modela, modelb, domain, times, var, percentile=95):
    grida = get_result(modela, domain, times)
    gridb = get_result(modelb, domain, times)
    return minmax([grida[f"{var}_{t}"] - gridb[f"{var}_{t}"] for t in times], percentile=percentile)
    
def compute_ranges(modelname: str, times: List[int], percentile=95):
    ranges = dict()
    ranges["sas"] = get_range(modelname, "sas", times, "c", percentile=percentile)
    ranges["arteries"] = get_range(modelname, "arteries", times, "c", percentile=percentile)
    ranges["veines"] = get_range(modelname, "venes", times, "c", percentile=percentile)

    return ranges

def compute_diff_ranges(modela:str, modelb:str, times: List[int],percentile=95):
    ranges = dict()
    ranges["sas"] = get_diff_range(modela, modelb, "sas", times, "c_sas", percentile=percentile)
    ranges["arteries"] = get_diff_range(modela, modelb, "arteries", times, "c_artery", percentile=percentile)
    ranges["veines"] = get_diff_range(modela, modelb, "venes", times, "c_vein", percentile=percentile)

    return ranges



def time_str(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return f"{h:02d}:{m:02d}"

def read_config(configfile):
    with open(configfile) as conf_file:
        config = yaml.load(conf_file, Loader=yaml.UnsafeLoader)
    return config


def clip_plot(csf, par, networks, filename, title, clim, cmap, cbar_title):
    import pyvista as pv
    csf_clipped = csf.clip(normal="y")
    par_clipped = par.clip(normal="y")
    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(csf_clipped, cmap=cmap, clim=clim,
                scalar_bar_args=dict(title=cbar_title, vertical=False,
                                    height=0.1, width=0.6, position_x=0.2,
                                    position_y=-0.0, title_font_size=36,
                                    label_font_size=32))
    pl.add_mesh(par_clipped, cmap=cmap, clim=clim, show_scalar_bar=False)
    for netw in networks:
        pl.add_mesh(netw, cmap=cmap, clim=clim, show_scalar_bar=False,
                    render_lines_as_tubes=True,line_width=5)

    pl.camera_position = 'zx'
    pl.camera.roll += 90
    pl.camera.zoom(1.3)
    pl.add_title(title, font_size=12)
    return pl.screenshot(filename, transparent_background=True, return_img=True)

def isosurf_plot(sas, networks, filename, title, clim, cmap, cbar_title):
    import pyvista as pv
    n = 5
    if clim is not None:
        contours = sas.contour(np.linspace(clim[0] + clim[1]/n, clim[1], n))
    else:
        contours = sas.contour(n)
        
    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(sas.outline(), color="white")
    pl.add_mesh(contours, opacity=0.75, clim=clim, cmap=cmap, show_scalar_bar=False)
    #pl.add_mesh(sas.extract_surface(), opacity=0.3, clim=clim, cmap=cmap,
    #            scalar_bar_args=dict(title=cbar_title, vertical=False,
    #                                height=0.1, width=0.6, position_x=0.2,
    #                                position_y=0.0, title_font_size=36,
    #                                label_font_size=32))
    for netw in networks:
        if clim is None: clim = (0, contours.active_scalars.max())
        netthres = netw.threshold(clim[1] / n)
        pl.add_mesh(netthres, cmap=cmap, clim=clim, show_scalar_bar=False,
                    render_lines_as_tubes=True,line_width=5)

    pl.camera_position = 'yz'
    #pl.camera.roll += 90
    pl.camera.zoom(1.3)
    pl.add_title(title, font_size=12)
    return pl.screenshot(filename, transparent_background=False, return_img=True)

def detail_plot(sas, networks, objects, filename, center, normal, clim, cmap, cbar_title):
    import pyvista as pv
    slice = sas.slice(origin=center, normal=normal)
    pl = pv.Plotter(off_screen=True)
    cmax = slice.active_scalars.max()
    for obj in objects:
        if "color" not in obj.array_names:
            cmax = max(cmax, obj.active_scalars.max()) 
    if clim is None: clim = (0, cmax)
    pl.add_mesh(slice, cmap=cmap, clim=clim, show_scalar_bar=False,
                scalar_bar_args=dict(title=cbar_title, vertical=False,
                                    height=0.1, width=0.6, position_x=0.2,
                                    position_y=0.05, title_font_size=36,
                                    label_font_size=32, color="grey"))
    for netw in networks:
        pl.add_mesh(netw, cmap=cmap, clim=clim, show_scalar_bar=False,
                    render_lines_as_tubes=True, line_width=20)
    for obj in objects:
        obj_args = dict(line_width=20, show_scalar_bar=False, cmap=cmap,)
        obj_args.update(dict(obj.field_data))
        pl.add_mesh(obj, clim=clim, **obj_args)

    #pl.add_legend_scale()
    #pl.add_axes()
    #pl.enable_parallel_projection()
    pl.camera_position = [center + normal, center, (0,0,1)]

    return slice.active_scalars.max(), pl.screenshot(filename, transparent_background=True, return_img=True)
