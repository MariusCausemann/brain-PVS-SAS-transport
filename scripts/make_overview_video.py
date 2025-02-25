# Make sure that dependencies are in place first:
#
# $ conda activate pvs_transport_env
#

#from plotting_utils import get_result, read_config
import pyvista as pv
import typer
from pathlib import Path
# import numpy as np
# from tqdm import tqdm
# from plot_csf_flow import from_k3d
# import k3d.colormaps.paraview_color_maps as pcm
# from extract_vessels import get_tubes
# import seaborn as sns

def make_overview_movie(datadir: str):

    # Read the images as NIFTI
    # Arteries
    #filename = Path(datadir, "data", "pcbi.1007073.s007.nii.gz")

    pl = pv.Plotter()
    pl.set_background('black')
    
    # Segmented T1
    filename = Path(datadir, "data", "T1_synthseg_robust.nii.gz")
    print("Reading image from ", filename)
    reader = pv.get_reader(filename) 
    img = reader.read()

    slices = img.slice_orthogonal()
    
    pl.add_mesh(slices, cmap="gray")

    pl.show()
    pl.screenshot("foo.png")
    
    print("Have plotted image")
    
if __name__ == "__main__":
    typer.run(make_overview_movie)
