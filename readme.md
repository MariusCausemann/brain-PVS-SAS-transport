
#### Installation and dependencies

To install all dependencies with mamba:

`mamba env create -f environment.yml`

This will create the environment 'pvs_transport_env'. Activate with:

`mamba activate pvs_transport_env`

To use conda instead of mamba, replace mamba with conda in the above.

Sidenote: to remove (for some reason or the other) your existing mamba environent, do

`mamba env remove -n pvs_transport_env`


#### Brain and vasculature imaging data 

Image data were downloaded from Hodneland et al, PLOS Comp. Bio, 2023 (https://doi.org/10.1371/journal.pcbi.1007073) and placed under data/. These data include: 

* Arteries: data/pcbi.1007073.s007.nii.gz
* Veins: data/pcbi.1007073.s008.nii.gz

In addition, we received additional T1 images from Erlend Hodneland, see data/T1_synthseg.nii.gz.

To have a quick look at the data, do e.g. the following (in Python) to plot a slice of the T1 data (NB: untested since data update).
```
import nibabel
import matplotlib.pyplot as plt
data = nibabel.load("data/T1_synthseg.nii.gz")
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

For networks/, there is the original, a smoothened and a tube representation of the arterial and venous network separately in .vtk format. The ones we use further are arteries_smooth.vtk and veins_smooth.vtk. 

For surfaces/, there is the white matter (pial) surface (wm.ply), the white-gray matter interface (gm.ply) and parenchyma.ply, all in .ply format. 

For volmesh/, there is the generated volumetric mesh of the parenchyma, including both white and gray matter (mesh.xdmf+h5).


#### Steps to run CSF flow on ex3
* login on ex3
* clone git repo
* install environment with conda/mamba - `mamba env create -f environment.yml`
* activate environment - `conda activate pvs_transport_env`
* run scrpt with slurm - `srun python scripts/sas_flow.py`

