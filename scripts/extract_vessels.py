import nibabel
import pyvista as pv
import numpy as np
from skimage.filters import frangi, sato
import skimage.morphology as skim
from skimage.measure import label
import kimimaro
import os


def np2pv(arr, resolution, origin):
    grid = pv.ImageData(dimensions=arr.shape + np.array((1, 1, 1)),
                        spacing=resolution, origin=origin)
    grid[f"data"] = arr.flatten(order="F")
    return grid

def skel_to_mesh(skel):
    n_segs = skel.edges.shape[0]
    cells = np.concatenate([np.ones((n_segs,1), dtype=int)*2, skel.edges], axis=1)
    celltypes = np.full(n_segs, fill_value=pv.CellType.LINE, dtype=np.uint8)
    mesh = pv.UnstructuredGrid(cells.ravel(), celltypes, skel.vertices)
    mesh["radius"] = skel.radius
    return mesh

def skeletonize(img):
    skels = kimimaro.skeletonize(img, anisotropy=(1, 1, 1),
                             teasar_params={"scale": scale, "const": const,})
    joined_skels = kimimaro.join_close_components(skels.values(), radius=50)
    ds_skel = joined_skels.downsample(2)
    return skel_to_mesh(ds_skel), ds_skel 

scale = 1.0
const = 1

os.makedirs("../mesh/networks", exist_ok=True)
# get the network centerlines of the venous network 
data = nibabel.load("../data/pcbi.1007073.s007.nii")
img = data.get_fdata().astype(int)
arterial_network, arteries = skeletonize(img)
with open('../mesh/arteries.swc','w') as f:
     f.write(arteries.to_swc())

arterial_network.save("../mesh/networks/arteries.vtk")

data = nibabel.load("../data/pcbi.1007073.s008.nii")
img = data.get_fdata().astype(int)
venous_network, veins = skeletonize(img)
with open('../mesh/veins.swc','w') as f:
     f.write(veins.to_swc())

venous_network.save("../mesh/networks/venes.vtk")


