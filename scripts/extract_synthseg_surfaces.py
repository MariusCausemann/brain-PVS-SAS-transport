import nibabel 
import pyvista as pv
import numpy as np
import skimage.morphology as skim
import os
import skimage
import scipy.ndimage as ndi
from scipy.ndimage import distance_transform_edt
from plotting_utils import read_config
from pathlib import Path
import typer

SYNTHSEG_V3 = 14
SYNTHSEG_V4 = 15
SYNTHSEG_LV = [4, 43]
SYNTHSEG_LV_INF = [5, 44]
SYNTHSEG_VENTRICLE = [SYNTHSEG_V3, SYNTHSEG_V4] +  SYNTHSEG_LV + SYNTHSEG_LV_INF
SYNTHSEG_CSF = [24]
SYNTHSEG_WM_LEFT, SYNTHSEG_WM_RIGHT = [2, 41]
SYNTHSEG_GM_LEFT, SYNTHSEG_GM_RIGHT,  = [3, 42]
SYNTHSEG_GM_CEREBELLUM_LEFT, SYNTHSEG_GM_CEREBELLUM_RIGHT = [8, 47]
SYTNSEG_CEREBRUM = [SYNTHSEG_WM_LEFT, SYNTHSEG_WM_RIGHT, SYNTHSEG_GM_LEFT, SYNTHSEG_GM_RIGHT]
SYNTHSEG_BRAIN_STEM = 16
mm2m = 1e-3

def get_closest_point(a, b):
    dist = ndi.distance_transform_edt(a==False)
    dist[b==False] = np.inf
    minidx = np.unravel_index(np.argmin(dist), a.shape)
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

def connect_by_line(m1, m2, footprint=skim.ball(1)):
    # compute connection between V3 and V4:
    pointa = get_closest_point(m1, m2)
    pointb = get_closest_point(m2, m1)

    # add a line between the shortest points to connect V3 and V4
    line = np.array(skimage.draw.line_nd(pointa, pointb, endpoint=True))
    conn = np.zeros(m1.shape, dtype=np.uint8)
    i,j,k = line
    conn[i,j,k] = 1
    return skim.binary_dilation(conn, footprint=footprint)


def seperate_labels(img, l1, l2, dist, newlabel=SYNTHSEG_CSF):
    m1 = skim.binary_dilation(img == l1, skim.ball(dist))
    m2 = skim.binary_dilation(img == l2, skim.ball(dist))
    img[np.logical_and(m1,m2)] = newlabel

def get_distance(img, labels):
    return distance_transform_edt(np.isin(img, labels)==False)

def create_tentorium(img):
    seperate_labels(img, SYNTHSEG_GM_CEREBELLUM_LEFT, SYNTHSEG_GM_LEFT, 3)
    seperate_labels(img, SYNTHSEG_GM_CEREBELLUM_RIGHT, SYNTHSEG_GM_RIGHT, 3)

    cerebellum_dist = get_distance(img, [SYNTHSEG_GM_CEREBELLUM_LEFT, 
                                         SYNTHSEG_GM_CEREBELLUM_RIGHT])
    cerebrum_dist = get_distance(img, [SYNTHSEG_GM_LEFT,SYNTHSEG_GM_RIGHT])
    brainstem_dist = get_distance(img, [SYNTHSEG_BRAIN_STEM])

    cerebellum_csf = np.logical_and(cerebellum_dist < cerebrum_dist, img==SYNTHSEG_CSF)
    cerebrum_csf = np.logical_and(cerebellum_dist > cerebrum_dist,  img==SYNTHSEG_CSF)

    tentorium_mask = np.logical_and(skim.binary_dilation(cerebellum_csf, skim.ball(1)), 
                                    skim.binary_dilation(cerebrum_csf, skim.ball(1)) )

    tentorium_mask[brainstem_dist < 8] = False
    tentorium_mask = skim.remove_small_objects(tentorium_mask, 1e3)
    return tentorium_mask


def extract_all_surfaces(configfile):

    meshname = Path(configfile).stem
    outdir = f"mesh/{meshname}/surfaces"
    os.makedirs(outdir, exist_ok=True)

    config = read_config(configfile)
    #load white matter data 
    pad_width = 5
    seg_rob = nibabel.load("data/T1_synthseg_robust.nii.gz")
    seg = nibabel.load("data/T1_synthseg.nii.gz")
    img_rob = np.pad(seg_rob.get_fdata(), pad_width=pad_width)
    img = np.pad(seg.get_fdata(), pad_width=pad_width)
    img[img_rob==0] = 0
    img[img_rob==SYNTHSEG_BRAIN_STEM] = SYNTHSEG_BRAIN_STEM
    resolution = np.array(seg.header["pixdim"][1:4])
    origin = - np.array(resolution) * pad_width

    hemisphere_min_distance = config.get("hemisphere_min_distance", 0)
    if hemisphere_min_distance:
        seperate_labels(img, SYNTHSEG_GM_LEFT, SYNTHSEG_GM_RIGHT,
                         hemisphere_min_distance)


    # generate skull surface - everything but background
    skull_mask = skim.remove_small_holes(img > 0, 1e3)
    par_mask = np.isin(img, SYNTHSEG_CSF + SYNTHSEG_VENTRICLE + [0]) == False
    skull_mask = np.logical_or(skull_mask, skim.binary_dilation(par_mask))
    skull_mask = skim.binary_dilation(skull_mask, skim.ball(1))

    if config.get("create_tentorium", False):
        tentorium_mask = create_tentorium(img)
        img[tentorium_mask] = 0
        skull_mask[tentorium_mask] = 0

    # adjust the bottom surface area - too small leads to to large pressure drop
    skull_mask = skim.binary_dilation(binary_smoothing(skim.binary_erosion(skull_mask)))
    skull_mask[:,:, :6] = np.logical_or(skull_mask,
        skim.binary_dilation(skull_mask, skim.ball(2)))[:,:, :6]
    skull_surf = extract_surface(skull_mask, resolution=resolution, origin=origin)
    skull_surf = skull_surf.smooth_taubin(n_iter=50, pass_band=0.1)
    #skull_surf.points = nibabel.affines.apply_affine(seg.affine, skull_surf.points)
    skull_surf = skull_surf.clip_closed_surface(normal=(0,0,1), origin=(0,0, 1))
    skull_surf.compute_normals(inplace=True, flip_normals=False)
    pv.save_meshio(f"{outdir}/skull.ply", skull_surf.scale(mm2m))

    # generate parenchyma surface- everything but CSF space
    par_mask = np.isin(img, SYNTHSEG_CSF + SYNTHSEG_VENTRICLE + [0]) == False
    par_mask = skim.remove_small_objects(par_mask, 1e3)
    par_mask = skim.remove_small_holes(par_mask, 1e3)
    par_surf = extract_surface(par_mask, resolution=resolution, origin=origin)
    par_surf = par_surf.smooth_taubin(n_iter=10, pass_band=0.1)
    #par_surf.points = nibabel.affines.apply_affine(seg.affine, par_surf.points)
    par_surf = par_surf.clip_closed_surface(normal=(0,0,1), origin=(0,0, 1))
    par_surf.compute_normals(inplace=True, flip_normals=False)
    pv.save_meshio(f"{outdir}/parenchyma.ply", par_surf.scale(mm2m))

    # generate cerebrum surface
    cerebrum_mask = np.isin(img, [SYTNSEG_CEREBRUM])
    cerebrum_mask = skim.remove_small_objects(cerebrum_mask, 1e3)
    cerebrum_mask = skim.remove_small_holes(cerebrum_mask, 1e3)
    cerebrum_surf = extract_surface(cerebrum_mask, resolution=resolution, origin=origin)
    cerebrum_surf = cerebrum_surf.smooth_taubin(n_iter=10, pass_band=0.1)
    cerebrum_surf.compute_normals(inplace=True, flip_normals=False)
    pv.save_meshio(f"{outdir}/cerebrum.ply", cerebrum_surf.scale(mm2m))

    # make inferior lateral ventricle horns

    for LVINFID, LVID in zip(SYNTHSEG_LV_INF, SYNTHSEG_LV):
        mask = img == LVINFID
        hull = skim.convex_hull_image(mask)
        mask = binary_smoothing(skim.binary_erosion(hull, skim.ball(1)) + 
                                skim.binary_dilation(mask, skim.ball(1)),
                                footprint=skim.ball(2))
        img[mask] = LVID
        surf = extract_surface(mask, resolution=resolution, origin=origin)
        surf = surf.smooth_taubin(n_iter=20, pass_band=0.05)
        pv.save_meshio(f"{outdir}/inf_{LVINFID}.ply", surf.scale(mm2m))

    # compute connection between V3 and V4:
    conn = connect_by_line(img == SYNTHSEG_V3, img==SYNTHSEG_V4, footprint=skim.ball(1.8))

    # compute ventricle surface
    ventricle_mask = np.isin(img, SYNTHSEG_VENTRICLE) + conn
    #ventricle_mask = skim.binary_dilation(ventricle_mask, footprint=skim.ball(1))
    ventricle_surf = extract_surface(ventricle_mask, resolution=resolution, origin=origin)
    ventricle_surf = ventricle_surf.smooth_taubin(n_iter=20, pass_band=0.05)
    ventricle_surf.compute_normals(inplace=True, flip_normals=False)
    pv.save_meshio(f"{outdir}/ventricles.ply", ventricle_surf.scale(mm2m))

    # compute lateral ventricle surface
    LV_mask = np.isin(img, SYNTHSEG_LV)
    for LVID in SYNTHSEG_LV:
        LV_mask += connect_by_line(img==LVID, img==SYNTHSEG_V3,footprint=skim.ball(2.5))
    #LV_mask =binary_smoothing(LV_mask, footprint=skim.ball(1))
    LV_mask = skimage.filters.gaussian(LV_mask)
    LV_surf = extract_surface(LV_mask, resolution=resolution, origin=origin)
    LV_surf = LV_surf.smooth_taubin(n_iter=20, pass_band=0.05)
    LV_surf.compute_normals(inplace=True, flip_normals=False)
    pv.save_meshio(f"{outdir}/LV.ply", LV_surf.scale(mm2m))

    # compute V3 and V4 surface
    V34_mask = np.isin(img, [SYNTHSEG_V3, SYNTHSEG_V4]) + conn
    V34_mask = skim.remove_small_objects(V34_mask, 1e3)
    V34_mask = skim.binary_dilation(V34_mask, footprint=skim.ball(1))
    V34_surf = extract_surface(V34_mask, resolution=resolution, origin=origin)
    V34_surf = V34_surf.smooth_taubin(n_iter=20, pass_band=0.05)
    V34_surf.compute_normals(inplace=True, flip_normals=False)
    pv.save_meshio(f"{outdir}/V34.ply", V34_surf.scale(mm2m))

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
    pv.save_meshio(f"{outdir}/parenchyma_incl_ventr.ply", par_surf.scale(mm2m))


if __name__ == "__main__":
    typer.run(extract_all_surfaces)