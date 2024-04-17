import nibabel 
import pyvista as pv
import numpy as np
import skimage.morphology as skim
import os
import typer
from plotting_utils import read_config
from pathlib import Path

def extract_surface(img):
    # img should be a binary 3D np.array
    resolution = (1,1,1)
    grid = pv.ImageData(dimensions=img.shape, spacing=resolution, origin=(0, 0, 0))
    mesh = grid.contour([0.5], img.flatten(order="F"), method="marching_cubes")
    surf = mesh.extract_geometry()
    return surf

def binary_smoothing(img, footprint=None):
    openend = skim.binary_opening(img, footprint=footprint)
    return skim.binary_closing(openend, footprint=footprint)

def extract_surfaces(configfile : str):

    config = read_config(configfile)
    meshname = Path(configfile).stem

    os.makedirs(f"mesh/{meshname}/surfaces", exist_ok=True)
    #load white matter data 
    wmdata = nibabel.load("data/pcbi.1007073.s005.nii.gz")
    wmimg = wmdata.get_fdata() 
    surf = extract_surface(wmimg)
    surf = surf.smooth_taubin(n_iter=20, pass_band=0.05)
    surf.save(f"mesh/{meshname}/surfaces/wm.ply")

    #load grey matter data 
    gmdata = nibabel.load("data/pcbi.1007073.s006.nii.gz")
    gmimg = gmdata.get_fdata() 
    surf_grey = extract_surface(gmimg)
    smooth_taubin_grey = surf_grey.smooth_taubin(n_iter=20, pass_band=0.05)
    smooth_taubin_grey.save(f"mesh/{meshname}/surfaces/gm.ply")

    # create parenchmya surface
    wmimg = np.pad(wmimg, ((0,0), (0,0), (0,8)))
    img = wmimg + gmimg
    for i in range(config["parenchyma_smooth_iterations"]):
        img = binary_smoothing(img, skim.ball(config["parenchyma_smooth_radius"]))
    surf_grey = extract_surface(img)
    smooth_taubin_parenchyma = surf_grey.smooth_taubin(n_iter=20, pass_band=0.05)
    smooth_taubin_parenchyma.save(f"mesh/{meshname}/surfaces/parenchyma.ply")

    ## creating a skull 
    ball = skim.ball(config["skull_dilate_radius"]) 
    img = np.pad(img, (0, 20))
    for i in range(config["skull_dilate_iterations"]):
        img = skim.binary_dilation(img, footprint=ball)
    img = skim.remove_small_holes(img, 1e6)
    for i in range(config["skull_smooth_iterations"]):
        img = binary_smoothing(img, skim.ball(config["skull_smooth_radius"]))
    for i in range(config["skull_erode_iterations"]):
        img = skim.erosion(img, skim.ball(config["skull_erode_radius"]))
    surf_dilated = extract_surface(img) 

    surf_dilated  = surf_dilated.smooth_taubin(n_iter=10, pass_band=0.05)
    surf_dilated.save(f"mesh/{meshname}/surfaces/skull.ply")


if __name__ == "__main__":
    typer.run(extract_surfaces)
