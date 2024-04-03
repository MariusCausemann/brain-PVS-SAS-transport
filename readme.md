
#### Installation and dependencies

To install dependencies with mamba (or conda):

`mamba env  create -f environment.yml`

#### Image data 

Image data were downloaded from Hodneland et al, PLOS Comp. Bio, 2023 (https://doi.org/10.1371/journal.pcbi.1007073) and placed under data/. These data include: 

* White matter representation (mask): data/pcbi.1007073.s005.nii.gz
* Gray matter representation (mask): data/pcbi.1007073.s006.nii.gz  
* Arteries: data/pcbi.1007073.s007.nii.gz
* Veins: data/pcbi.1007073.s008.nii.gz

To have a quick look at the data, do e.g. the following (in Python) to plot a slice of the white matter data.
```
import nibabel
import matplotlib.pyplot as plt
data = nibabel.load("data/pcbi.1007073.s005.nii.gz")
data.shape
plt.imshow(data.get_fdata()[:, 100, :])
plt.show()
```

#### Generating surfaces, meshes and networks from the image data

To generate surfaces, extract vessel networks and generate meshes, run:

```
cd scripts
python3 extract_surfaces.py
python3 extract_vessels.py
python3 generate_mesh.py
```

These scripts will generate objects placed under mesh/ including

* networks/
* surfaces/
* volmesh/

For networks/, there is the original, a smoothened and a tube representation of the arterial and venous network separately in .vtk format.

For surfaces/, there is the white matter (pial) surface (wm.ply), the white-gray matter interface (gm.ply) and parenchyma.ply, all in .ply format. 

For volmesh/, there is the generated volumetric mesh of the parenchyma, including both white and gray matter (mesh.xdmf+h5).
