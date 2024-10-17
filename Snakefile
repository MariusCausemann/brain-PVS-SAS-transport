import numpy as np
from scripts.plotting_utils import read_config
from collections import defaultdict

num_mumps_threads = 16
mesh_refine_models = ["modelALowRes", "modelA", "modelAHighRes"]
time_refine_models = ["modelA", "modelAsmalldt", "modelAlargedt"]

models =  mesh_refine_models + time_refine_models
times = list(np.array([1, 6, 12, 24])*60*60)
conctimes =  list((np.array([0, 1/3, 2/3, 1, 4/3, 5/3, 2, 7/3, 8/3, 3, 4, 5, 6 ,9, 12, 15, 18, 21, 24])*60*60).astype(int))

cmax = {"detail":defaultdict(lambda: 2, {"modelA_modelA4":2}),
        "overview":defaultdict(lambda: 2, {"modelA_modelA4":2}),
        "isosurf":defaultdict(lambda: 5, {"modelA_modelA4":8}),
        }

diffmax = {"detail":defaultdict(lambda: 0.1,{"modelA_modelA4":0.2}),
           "overview":defaultdict(lambda: 1,{"modelA_modelA4":1}),
           "isosurf":defaultdict(lambda: 1,{"modelA_modelA4":1}),          
}
types = ["overview", "isosurf"]

# check for headless server here? https://github.com/pyvista/pyvista/blob/f872ffebc7f56b1ff03541e368935c3304fc5e33/pyvista/plotting/tools.py#L54
# or try vtk osmesa build...

def getconfig(model, key):
    return read_config(f"configfiles/{model}.yml").get(key, [])

rule all:
    input:
        #"plots/comparisons/modelA_modelALowRes/modelA_modelALowRes_overview.png",
        "plots/comparisons/modelA_modelAHighRes/modelA_modelAHighRes_overview.png",
        "plots/comparisons/modelA_modelAHighRes/modelA_modelAHighRes_isosurf.png",
        #"plots/comparisons/modelA_modelA3/modelA_modelA3_overview.png",
        #"plots/comparisons/modelA_modelA4/modelA_modelA4_overview.png",
        #"plots/comparisons/modelA_modelA3/modelA_modelA3_overview.png",
        #"plots/comparisons/modelB_modelE/modelB_modelE_isosurf.png",
        #"plots/comparisons/modelA_modelB/modelA_modelB_detail.png",
        expand("plots/{modelname}/{modelname}_tracer_vessel_dist.png", modelname=models),
        expand("plots/{modelname}/{modelname}_total_conc.png", modelname=models),
        expand("plots/{modelname}/{modelname}_{tp}_{t}.png", modelname=models, t=times, tp=types)


rule runSimuation:
    conda:"environment.yml"
    input:
        volmesh=lambda wildcards: getconfig(wildcards.modelname, "mesh"),
        artmesh="mesh/networks/arteries_smooth.vtk",
        venmesh="mesh/networks/venes_smooth.vtk",
        config="configfiles/{modelname}.yml",
        csf_velocity_file=lambda wildcards: getconfig(wildcards.modelname, "csf_velocity_file"),
        arterial_velocity_file=lambda wildcards: getconfig(wildcards.modelname, "arterial_velocity_file"),
        csf_dispersion_file=lambda wildcards: getconfig(wildcards.modelname, "csf_dispersion_file"),
    output:
        sas="results/{modelname}/{modelname}_sas.xdmf",
        art="results/{modelname}/{modelname}_artery.xdmf",
        ven="results/{modelname}/{modelname}_vein.xdmf",
    threads: num_mumps_threads
    resources:
        ncpuspertask=num_mumps_threads
    shell:
        "export OMP_NUM_THREADS={threads} && python scripts/time_dependent.py {input.config}"

rule computeDispersionField:
    conda:"environment.yml"
    input:
        csf_pressure_file="results/csf_flow/{csf_flow_model}/flow.hdf",
        bg="mesh/standard/standard.xdmf"
    output:
        csf_dispersion_file="results/csf_flow/{csf_flow_model}/R.xdmf",
        csf_dispersion_plot="results/csf_flow/{csf_flow_model}/R.png",
    shell:
        "xvfb-run -a python scripts/compute_dispersion_field.py {input.csf_pressure_file} {output.csf_dispersion_file}"

rule computeFlowField:
    conda:"environment.yml"
    input:
        artmesh="mesh/networks/arteries_smooth.vtk",
    output:
        flowfield="results/pvs_flow/pvs_flow.xdmf",
    shell:
        "python scripts/pvs_flow.py"

rule generatePlot:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.xdmf",
        art="results/{modelname}/{modelname}_artery.xdmf",
        ven="results/{modelname}/{modelname}_vein.xdmf",
    output:
        plot="plots/{modelname}/{modelname}_{type}_{time,[0-9]*}.png"
    shell:
        "xvfb-run -a python scripts/generate_plot.py {wildcards.modelname} {wildcards.time} {wildcards.type} --filename {output.plot}"

rule generateDiffPlot:
    conda:"environment.yml"
    input:
        sas1="results/{model1}/{model1}_sas.xdmf",
        art1="results/{model1}/{model1}_artery.xdmf",
        ven1="results/{model1}/{model1}_vein.xdmf",
        sas2="results/{model2}/{model2}_sas.xdmf",
        art2="results/{model2}/{model2}_artery.xdmf",
        ven2="results/{model2}/{model2}_vein.xdmf",
    output:
        plot="plots/comparisons/{model1}_{model2}/{model1}_{model2}_diff_{type}_{time}.png"
    shell:
        "xvfb-run -a python scripts/create_diff_plot.py {wildcards.model1} {wildcards.model2} {wildcards.time} {wildcards.type} --filename {output.plot}"

rule getConcentrationRange:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.xdmf",
        art="results/{modelname}/{modelname}_artery.xdmf",
        ven="results/{modelname}/{modelname}_vein.xdmf",
    output:
        ranges="plots/{modelname}/{modelname}_ranges.yml"
    params:
        times=times
    shell:
        "python scripts/compute_ranges.py {wildcards.modelname} {params.times}"


rule compareModels:
    conda:"environment.yml"
    input:
        sas1="results/{model1}/{model1}_sas.xdmf",
        art1="results/{model1}/{model1}_artery.xdmf",
        ven1="results/{model1}/{model1}_vein.xdmf",
        sas2="results/{model2}/{model2}_sas.xdmf",
        art2="results/{model2}/{model2}_artery.xdmf",
        ven2="results/{model2}/{model2}_vein.xdmf"
    output:
        plot="plots/comparisons/{model1}_{model2}/{model1}_{model2}_{type}.png"
    params:
        cmax= lambda wildcards: cmax[wildcards.type][f"{wildcards.model1}_{wildcards.model2}"],
        diffmax= lambda wildcards: diffmax[wildcards.type][f"{wildcards.model1}_{wildcards.model2}"]
    shell:
        """
        xvfb-run -a python scripts/compare_models.py {wildcards.model1} {wildcards.model2} {wildcards.type} {params.cmax} {params.diffmax} {times} &&
        xvfb-run -a python scripts/compare_models_horizontal.py {wildcards.model1} {wildcards.model2} {wildcards.type} {params.cmax} {times}
        """

rule analyzeTracerDist:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.xdmf",
        art="results/{modelname}/{modelname}_artery.xdmf",
        ven="results/{modelname}/{modelname}_vein.xdmf",
    output:
        plot="plots/{modelname}/{modelname}_tracer_vessel_dist.png"
    shell:
        "xvfb-run -a python scripts/analyze_tracer_dist.py {wildcards.modelname} {times}"

rule totalTracer:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.xdmf",
        art="results/{modelname}/{modelname}_artery.xdmf",
        ven="results/{modelname}/{modelname}_vein.xdmf",
    output:
        plot="plots/{modelname}/{modelname}_total_conc.png"
    shell:
        "xvfb-run -a python scripts/mean_concentrations.py {wildcards.modelname} {conctimes}"

rule segmentT1:
    input:
        "data/T1.nii.gz"
    output:
        "data/T1_synthseg.nii.gz"
    shell:
        "scripts/synthseg.sh"


rule generateSurfaces:
    conda:"environment.yml"
    input:
        "data/T1_synthseg.nii.gz"
    output:
        [f"mesh/T1/surfaces/{n}.ply" for n in ["LV", "parenchyma", "skull", "V34"]]
    shell:
        """
        python scripts/extract_synthseg_surfaces.py
        """

rule generateMesh:
    conda:"environment.yml"
    input:
        [f"mesh/T1/surfaces/{n}.ply" for n in ["LV", "parenchyma", "skull", "V34"]],
        "configfiles/meshconfig/{meshname}.yml"
    output:
        "mesh/{meshname}/volmesh/mesh.xdmf",
    resources:
        ncpuspertask=4
    shell:
        """
        python scripts/generate_synthseg_mesh.py configfiles/meshconfig/{wildcards.meshname}.yml
        """

rule markAndRefineMesh:
    conda:"environment.yml"
    input:
        "mesh/{meshname}/volmesh/mesh.xdmf",
    output:
        "mesh/{meshname}/{meshname}.xdmf",
        "mesh/{meshname}/{meshname}.h5",
        "mesh/{meshname}/{meshname}_outer.xdmf",
        "mesh/{meshname}/{meshname}_outer.h5",
    shell:
        """
        python scripts/mark_and_refine_mesh.py configfiles/meshconfig/{wildcards.meshname}.yml
        """

rule computeSASFlow:
    conda:"environment.yml"
    input:
        volmesh=lambda wildcards: getconfig(wildcards.csf_flow_model, "mesh"),
        config="configfiles/{csf_flow_model}.yml",
    output:
        "results/csf_flow/{csf_flow_model}/csf_p.xdmf",
        "results/csf_flow/{csf_flow_model}/csf_p.h5",
        "results/csf_flow/{csf_flow_model}/csf_v.xdmf",
        "results/csf_flow/{csf_flow_model}/csf_v.h5",
        "results/csf_flow/{csf_flow_model}/flow.hdf",
    threads: 4
    resources:
        ncpuspertask=64
    shell:
        """
        export OMP_NUM_THREADS={threads} && \
        mpiexec -n 16 python scripts/sas_flow_Hdiv.py \
        configfiles/{wildcards.csf_flow_model}.yml && \
        xvfb-run -a python3 scripts/plot_csf_flow.py results/csf_flow/{wildcards.csf_flow_model}
        """


rule computePVSFlow:
    conda:"environment.yml"
    input:
        "results/csf_flow/{csf_flow_model}/flow.hdf",
    output:
        'results/pvs_flow_prod/{csf_flow_model}/pvs_flow.xdmf',
        'results/pvs_flow_prod/{csf_flow_model}/pvs_flow.h5'
    shell:
        "python scripts/pvs_flow_prod.py {wildcards.csf_flow_model}"


rule AverageSASFlow2PVS:
    conda:"environment.yml"
    input:
        "results/csf_flow/{csf_flow_model}/flow.hdf",
    output:
        'results/csf_flow/{csf_flow_model}/pvs_flow.xdmf',
        'results/csf_flow/{csf_flow_model}/pvs_flow.h5'
    shell:
        "python scripts/uhat_prod.py {wildcards.csf_flow_model}"





    







