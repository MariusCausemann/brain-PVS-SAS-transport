![Actions Status](https://github.com/mariuscausemann/brain-PVS-SAS-transport/actions/workflows/test_conda.yml/badge.svg)

#### Installation and dependencies

The workflow is managed with the workflow management system [Snakemake](https://snakemake.readthedocs.io/en/stable/). To install snakemake using conda, run:

`conda install -c conda-forge -c bioconda snakemake==8.14.0`

Next, you can test the workflow with a small scale example with:

`snakemake --conda-frontend conda --use-conda --cores 2 -p plots/test/test_total_conc.png --config meshing=False`

This will automatically install all required dependencies (from environment.yml), and run all required jobs on two cores. Since the meshing tool fTetWild requires compilation and can be hard to install, we disable the meshing part of the pipeline here and run the example on a pregenerated mesh. Expected run time for setting up the environment and running all steps is ~ 30min.
To test the whole pipline (including meshing), run:

`snakemake --conda-frontend conda --use-conda --cores 2 -p plots/test/test_total_conc.png --force-all`

To reproduce all results on N cores, run:

`snakemake --conda-frontend conda --use-conda --cores N`

Note that this requires significant computational resources (around ~ 500GB RAM, 64 cores minimum). Snakemake supports job submission systems like `slurm`. Adjust our profile in `ex3/config.yaml` to your needs, install the generic executor plugin with

`conda install -c conda-forge -c bioconda snakemake-executor-plugin-cluster-generic`

and start with:

`snakemake --profile ex3`

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

#### Generating meshes from the image data ####

The standard mesh resolution is available in the Git repository under
mesh/standard.


For networks/, there is the original, a smoothened and a tube representation of the arterial and venous network separately in .vtk format. The ones we use further are arteries_smooth.vtk and veins_smooth.vtk. 

For surfaces/, there is the white matter (pial) surface (wm.ply), the white-gray matter interface (gm.ply) and parenchyma.ply, all in .ply format. 

For volmesh/, there is the generated volumetric mesh of the parenchyma, including both white and gray matter (mesh.xdmf+h5).
