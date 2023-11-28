import pyvista as pv
from plotting_utils import get_result, read_config
from pykdtree.kdtree import KDTree
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from plotting_utils import time_str, set_plotting_defaults
import typer


def plot_concentration_distance(c_time_series, isodists, filename):
    set_plotting_defaults()
    plt.figure()
    for t, c in c_time_series.items():
        plt.plot(isodists, c, label=f"{time_str(t)} h", marker="o")
    plt.legend()
    plt.xlabel("vessel distance (m)")
    plt.ylabel("concentration")
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(filename)


def analyze_tracer_dist(modelname: str):
    domain = "sas"
    times = 3600*np.array([1, 6, 12, 18, 24])
    sas = get_result(modelname, domain, times)
    art = get_result(modelname, "arteries", 0)
    ven = get_result(modelname, "venes", 0)
    tree = KDTree(np.concatenate((art.points, ven.points), axis=0))
    d, idx = tree.query(sas.points)
    sas["distances"] = d
    isodists = np.linspace(0, 0.1, 10)
    isosurfs = [sas.contour(isosurfaces=[d], scalars="distances") for d in isodists]

    c_time_series = {}
    for t in times:
        c_time_series[t] = [isos[f"c_sas_{t}"].mean() for isos in isosurfs]
    
    filename = f"plots/{modelname}/{modelname}_tracer_vessel_dist.png"
    plot_concentration_distance(c_time_series, isodists, filename)

if __name__ == "__main__":
    typer.run(analyze_tracer_dist)

