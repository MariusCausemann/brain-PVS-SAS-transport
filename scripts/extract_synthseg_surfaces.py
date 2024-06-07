import nibabel 
import pyvista as pv
import numpy as np
import skimage.morphology as skim
import os
import skimage
import scipy.ndimage as ndi


SYNTHSEG_V3 = 14
SYNTHSEG_V4 = 15
SYNTHSEG_LV = [4, 43]
SYNTHSEG_VENTRICLE = [SYNTHSEG_V3, SYNTHSEG_V4] +  SYNTHSEG_LV
SYNTHSEG_CSF = [24]

mm2m = 1e-3

def get_closest_point(a, b):
    dist = ndi.distance_transform_edt(a==False)
    dist[b==False] = np.inf
    minidx = np.unravel_index(np.argmin(dist), img.shape)
    return minidx

def extract_surface(img, resolution=(1,1,1), origin=(0, 0, 0)):
    # img should be a binary 3D np.array
    grid = pv.ImageData(dimensions=img.shape, spacing=resolution, origin=origin)
    mesh = grid.contour([0.5], img.flatten(order="F"), method="marching_cubes")
    surf = mesh.extract_geometry()
    surf.clear_data()
    return surf

def binary_smoothing(img, footprint=skim.ball(1)):
    openend = skim.binary_opening(img, footprint=footprint)
    return skim.binary_closing(openend, footprint=footprint)

os.makedirs("mesh/T1/surfaces", exist_ok=True)
#load white matter data 
pad_width = 5
seg = nibabel.load("results/freesurfer/T1_synthseg.nii.gz")
img = np.pad(seg.get_fdata(), pad_width=pad_width)
resolution = np.array(seg.header["pixdim"][1:4])
origin = - np.array(resolution) * pad_width


# generate skull surface - everything but background
skull_mask = img > 0
skull_mask = skim.binary_dilation(skull_mask, skim.ball(1))
skull_surf = extract_surface(skull_mask, resolution=resolution, origin=origin)
skull_surf = skull_surf.smooth_taubin(n_iter=20, pass_band=0.05)
#skull_surf.points = nibabel.affines.apply_affine(seg.affine, skull_surf.points)
skull_surf = skull_surf.clip_closed_surface(normal=(0,0,1), origin=(0,0, 1))
skull_surf.compute_normals(inplace=True, flip_normals=False)
pv.save_meshio("mesh/T1/surfaces/skull.ply", skull_surf.scale(mm2m))

# generate parenchyma surface- everything but CSF space
par_mask = np.isin(img, SYNTHSEG_CSF + SYNTHSEG_VENTRICLE + [0]) == False
par_surf = extract_surface(par_mask, resolution=resolution, origin=origin)
par_surf = par_surf.smooth_taubin(n_iter=20, pass_band=0.05)
#par_surf.points = nibabel.affines.apply_affine(seg.affine, par_surf.points)
par_surf = par_surf.clip_closed_surface(normal=(0,0,1), origin=(0,0, 1))
par_surf.compute_normals(inplace=True, flip_normals=False)
pv.save_meshio("mesh/T1/surfaces/parenchyma.ply", par_surf.scale(mm2m))

# compute connection between V3 and V4:
V3 = np.isin(img, SYNTHSEG_V3)
V4 = np.isin(img, SYNTHSEG_V4)
pointa = get_closest_point(V3, V4)
pointb = get_closest_point(V4, V3)

# add a line between the shortest points to connect V3 and V4
line = np.array(skimage.draw.line_nd(pointa, pointb, endpoint=True))
conn = np.zeros(img.shape)
i,j,k = line
conn[i,j,k] = 1
conn = skim.binary_dilation(conn, footprint=skim.ball(1))

# compute ventricle surface
ventricle_mask = np.isin(img, SYNTHSEG_VENTRICLE) + conn
ventricle_mask = skim.binary_dilation(ventricle_mask, footprint=skim.ball(1))
ventricle_surf = extract_surface(ventricle_mask, resolution=resolution, origin=origin)
ventricle_surf = ventricle_surf.smooth_taubin(n_iter=20, pass_band=0.05)
ventricle_surf.compute_normals(inplace=True, flip_normals=False)
pv.save_meshio("mesh/T1/surfaces/ventricles.ply", ventricle_surf.scale(mm2m))

# compute lateral ventricle surface
LV_mask = np.isin(img, SYNTHSEG_LV)
LV_mask = skim.binary_dilation(LV_mask, footprint=skim.ball(1))
LV_mask = extract_surface(LV_mask, resolution=resolution, origin=origin)
LV_mask = LV_mask.smooth_taubin(n_iter=20, pass_band=0.05)
LV_mask.compute_normals(inplace=True, flip_normals=False)
pv.save_meshio("mesh/T1/surfaces/LV.ply", LV_mask.scale(mm2m))

# compute V3 and V4 surface
V34_mask = np.isin(img, [SYNTHSEG_V3, SYNTHSEG_V4]) + conn
V34_mask = skim.binary_dilation(V34_mask, footprint=skim.ball(1))
V34_mask = extract_surface(V34_mask, resolution=resolution, origin=origin)
V34_mask = V34_mask.smooth_taubin(n_iter=20, pass_band=0.05)
V34_mask.compute_normals(inplace=True, flip_normals=False)
pv.save_meshio("mesh/T1/surfaces/V34.ply", V34_mask.scale(mm2m))

# compute a layer of parenchymal tissue around the ventricles
# to get a watertight ventricular system and
# generate combined mask of ventricles and parenchyma

# first find lowest point of V4 to generate an outlet of the tissue sheet into
# the cisterna magna
ventr_idx = np.argwhere(ventricle_mask)
ind = np.argsort(ventr_idx[:, 2])
lowest_ventr_point = ventr_idx[ind][0]
ventr_outlet = np.zeros(img.shape)
ventr_outlet[tuple(lowest_ventr_point)] = 1
ventr_outlet = skim.binary_dilation(ventr_outlet, footprint=skim.ball(5))

# extent ventricles to create tissue sheet and substract cisterna magna outlet
ventricle_extended = skim.binary_dilation(ventricle_mask, footprint=skim.ball(3))
ventricle_extended = np.logical_and(ventricle_extended, ventr_outlet==False)

# finally, generate combined mask of parenchyma and ventricles 
par_ventr_mask = par_mask + ventricle_extended
par_surf = extract_surface(par_ventr_mask, resolution=resolution, origin=origin)
par_surf = par_surf.smooth_taubin(n_iter=20, pass_band=0.05)
#par_surf.points = nibabel.affines.apply_affine(seg.affine, par_surf.points)
par_surf = par_surf.clip_closed_surface(normal=(0,0,1), origin=(0,0, 1))
par_surf.compute_normals(inplace=True, flip_normals=False)
pv.save_meshio("mesh/T1/surfaces/parenchyma_incl_ventr.ply", par_surf.scale(mm2m))