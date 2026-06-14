![Actions Status](https://github.com/mariuscausemann/brain-PVS-SAS-transport/actions/workflows/test_conda.yml/badge.svg)

# brain-PVS-SAS-transport

This repository contains the computational framework for the paper:
**"In-silico solute transport via perivascular networks in the human intracranial space"** (Causemann et al., 2026).

This software provides a mixed-dimensional mathematical modelling framework to predict the spatiotemporal spreading of a solute across the subarachnoid space (SAS), ventricular system, and brain parenchyma, including surface perivascular spaces (PVSs). The simulation pipeline is fully automated using [Snakemake](https://snakemake.readthedocs.io/) and uses the [FEniCS](https://fenicsproject.org/) finite element framework for numerical approximations.

## 1. System Requirements

### Hardware Requirements

* **Minimal Demo:** A standard desktop computer with at least 2 CPU cores and 8 GB of RAM.

* **Full Reproduction:** A High-Performance Computing (HPC) environment is **required** to run the full three-dimensional Stokes flow and solute transport models.

  * *Stokes flow simulation:* ~16 MPI processes (64 threads) and ~300 GB RAM.

  * *Transport simulation:* ~16 threads and ~50 GB RAM.

### Software Dependencies & Operating Systems

* **Operating Systems:** Linux (tested on Ubuntu-latest via GitHub Actions). It should also run on macOS or Windows via WSL.

* **Core Dependencies:**

  * Conda / Miniconda / Mamba

  * Snakemake `==8.14.0` (with plugins `snakemake-storage-plugin-http` and `snakemake-executor-plugin-cluster-generic`)

* **Environment Dependencies:** All simulation-specific dependencies (e.g., Python 3.12, FEniCS, NumPy, YAML) are automatically resolved and installed via Conda environments defined in the repository (e.g., `environment.yml`, `mesh_environment.yml`) when running the Snakemake workflow. A list of the exact versions used can be found in `frozen_env.yml`.

## 2. Installation Guide

### Instructions

1. **Clone the repository:**

   ```
   git clone [https://github.com/MariusCausemann/brain-PVS-SAS-transport.git](https://github.com/MariusCausemann/brain-PVS-SAS-transport.git)
   cd brain-PVS-SAS-transport
   
   ```

2. **Install Snakemake via Conda/Mamba:**
   If you do not have Conda installed, please install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) first. Then, create an environment for Snakemake:

   ```
   conda create -n snakemake_env -c conda-forge -c bioconda python=3.12 git snakemake==8.14.0 snakemake-storage-plugin-http snakemake-executor-plugin-cluster-generic
   conda activate snakemake_env
   
   ```

### Typical Install Time

* Snakemake installation: ~2-5 minutes.

* Sub-environment creation (handled automatically by Snakemake upon first execution): ~5-10 minutes depending on internet speed.

## 3. Demo

To verify that the pipeline is functioning correctly, you can run a minimal demo. This demo executes a small subset of the pipeline to produce a basic tracer concentration plot, skipping the computationally expensive 3D meshing steps.

### Instructions to run on data

Ensure your `snakemake_env` is activated, then run:

```
snakemake --conda-frontend conda --use-conda --cores 2 -p plots/test/test_total_conc.png --config meshing=False

```

### Expected Output

The workflow will download/configure necessary dependencies via conda, run the basic numerical scripts, and generate an output plot located at:
`plots/test/test_total_conc.png`
Note that this is intended to test the setup, and will not produce physically meaningful simulation results.

### Expected Run Time

On a standard desktop computer, the minimal demo takes approximately **15-30 minutes** (including the time required for Conda to build the inner `environment.yml` for the first time). Subsequent runs will be faster.

## 4. Instructions for Use

### How to run the software on your data

The entire pipeline—from MRI segmentation extraction to mesh generation, Stokes flow calculation, mixed-dimensional transport simulation, and final plotting—is defined in the `Snakefile`.

You can configure the model runs by modifying the YAML files located in the `configfiles/` directory.

### Reproduction Instructions

To reproduce the quantitative results and figures presented in the manuscript, you must run the full workflow.
*Note: This will trigger the full FEniCS simulations on the high-resolution intracranial meshes and requires an HPC environment.*

1. **Dry-run (Optional):** See what rules will be executed without running them:

   ```
   snakemake -n
   
   ```

2. **Full Execution:** Run the entire pipeline utilizing your cluster/HPC resources (adjust the `--cores` flag based on your available hardware):

   ```
   snakemake --use-conda --cores 64
   
   ```

**Key Pipeline Steps included in the Workflow:**

* `segmentT1` / `generateSurfaces` / `generateMesh`: Extracts surfaces from T1-weighted MR images and generates the 3D computational mesh (fTetWild).

* `computeSASFlow` / `computeProdPVSFlow`: Computes steady and peristaltic flow fields in the Subarachnoid Space and Perivascular Spaces.

* `runSimuation`: Executes the mixed-dimensional 3D-1D time-dependent solute transport models.

* `generatePlot` / `compareModels` / `makeVideo`: Generates the resulting concentration visualizations, bar plots, and animations across 1h to 24h timelines.
