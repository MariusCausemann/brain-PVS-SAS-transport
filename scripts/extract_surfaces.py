import nibabel 
import pyvista as pv
import numpy as np
import skimage.morphology as skim
import os

def extract_surface(img):
    # img should be a binary 3D np.array
    resolution = (1,1,1)
    grid = pv.ImageData(dimensions=img.shape, spacing=resolution, origin=(0, 0, 0))
    mesh = grid.contour([0.5], img.flatten(order="F"), method="marching_cubes")
    surf = mesh.extract_geometry()
    return surf

def skel_to_mesh(skel):
    n_segs = skel.edges.shape[0]
    cells = np.concatenate([np.ones((n_segs,1), dtype=int)*2, skel.edges], axis=1)
    celltypes = np.full(n_segs, fill_value=pv.CellType.LINE, dtype=np.uint8)
    mesh = pv.UnstructuredGrid(cells.ravel(), celltypes, skel.vertices)
    mesh["radius"] = skel.radius
    return mesh

def binary_smoothing(img, footprint=None):
    openend = skim.binary_opening(img, footprint=footprint)
    return skim.binary_closing(openend, footprint=footprint)

os.makedirs("../mesh/surfaces", exist_ok=True)
#load white matter data 
wmdata = nibabel.load("../data/pcbi.1007073.s005.nii")
wmimg = wmdata.get_fdata() 
surf = extract_surface(wmimg)
surf = surf.smooth_taubin(n_iter=20, pass_band=0.05)
surf.save("../mesh/surfaces/wm.ply")

#load grey matter data 
gmdata = nibabel.load("../data/pcbi.1007073.s006.nii")
gmimg = gmdata.get_fdata() 
surf_grey = extract_surface(gmimg)
smooth_taubin_grey = surf_grey.smooth_taubin(n_iter=20, pass_band=0.05)
smooth_taubin_grey.save("../mesh/surfaces/gm.ply")

# create parenchmya surface
wmimg = np.pad(wmimg, ((0,0), (0,0), (0,8)))
img = wmimg + gmimg
surf_grey = extract_surface(img)
smooth_taubin_parenchyma = surf_grey.smooth_taubin(n_iter=20, pass_band=0.05)
smooth_taubin_parenchyma.save("../mesh/surfaces/parenchyma.ply")

## creating a skull 
ball = skim.ball(5) 
img = np.pad(img, (0, 20))
for i in range(3):
    img = skim.binary_dilation(img, footprint=ball)
img = skim.remove_small_holes(img, 1e6)
for i in range(3):
    img = binary_smoothing(img, skim.ball(5))
#img = skim.erosion(img, skim.ball(5))
surf_dilated = extract_surface(img) 

surf_dilated  = surf_dilated.smooth_taubin(n_iter=10, pass_band=0.05)
surf_dilated.save("../mesh/surfaces/skull.ply")