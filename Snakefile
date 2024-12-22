import numpy as np
from collections import defaultdict
import yaml

def read_config(configfile):
    with open(configfile) as conf_file:
        config = yaml.load(conf_file, Loader=yaml.UnsafeLoader)
    return config

num_mumps_threads = 16
mesh_refine_models = ["modelALowRes", "modelA", "modelAHighRes"]
time_refine_models = ["modelA", "modelAsmalldt", "modelAlargedt"]
model_variations = ["modelA" , "modelA-LowD","modelA-HighD",
                    "modelA-OnlyDispersion",
                    "modelA-strongVM","modelA-PVS-disp",
                    "modelB1-10", "modelB1-100", "modelB1-1000",
                    "modelB2-10", "modelB2-100", "modelB2-1000",
                    #"modelC", "modelA-NoResp", "modelA-NoDisp", "modelA-LowD",
                    "LargePVS", "LargePVSA"]

model_comparisons = [
                    #"modelA_modelB1-10","modelA_modelB1-100", "modelA_modelB1-1000",
                    #"modelA_modelB2-1", "modelA_modelB2-10", "modelA_modelB2-100",
                    # "modelA_modelB2-1000", "modelA_modelC",
                    # "modelA_modelA-NoResp",
                    "modelA_modelA-strongVM",
                    "modelA_modelA-PVS-disp"            
                    ]

models =  model_variations #mesh_refine_models + time_refine_models
times = list(np.array([1, 3, 6, 12, 24])*60*60)
conctimes =  list((np.array([0, 1/3, 2/3, 1, 4/3, 5/3, 2, 7/3, 8/3, 3, 4, 5, 6 ,9, 12, 15, 18, 21, 24])*60*60).astype(int))

cmax = {"detail":defaultdict(lambda: 2, {"modelA_modelA4":2}),
        "overview":defaultdict(lambda: 2, {"modelA_modelA4":2}),
        "isosurf":defaultdict(lambda: 5, {"modelA_modelA4":8}),
        "timesurf":defaultdict(lambda: 1, {"modelA_modelA4":1}),
        }

diffmax = {"detail":defaultdict(lambda: 0.1,{"modelA_modelA4":0.2}),
           "overview":defaultdict(lambda: 1,{"modelA_modelA4":1}),
           "isosurf":defaultdict(lambda: 1,{"modelA_modelA4":1}),  
           "timesurf":defaultdict(lambda: 1, {"modelA_modelA4":1}),        
}
types = ["overview", "isosurf", "timesurf"]

def getconfig(model, key):
    return read_config(f"configfiles/{model}.yml").get(key, [])

rule all:
    input:
        expand("plots/comparisons/{c}/{c}_{type}.png",c=model_comparisons, type=types),
        expand("plots/{modelname}/{modelname}_tracer_vessel_dist.png", modelname=models),
        expand("plots/{modelname}/{modelname}_total_conc.png", modelname=models),
        expand("plots/{modelname}/{modelname}_overview_1-6-12-24.png", modelname=models),
        expand("plots/{modelname}/{modelname}_overview_1-3-6-9-12-24.png", modelname=models),
        expand("plots/{modelname}/{modelname}_overview_4-6.png", modelname=models),
        expand("plots/{modelname}/{modelname}.mp4", modelname=models),
        expand("plots/{modelname}/{modelname}_1+2+3+6+9+12_all_0_details.png", modelname=models),
        "plots/pvs_flow_prod/sas_flow-arteries/",
        "plots/pvs_flow_peristaltic/vasomotion/",
        "plots/pvs_flow_peristaltic/cardiac_pvs_oscillation/",
        #expand("plots/{modelname}/{modelname}_{tp}_{t}.png", modelname=models, t=times, tp=types)
        expand("plots/{m}/{m}_{tstr}_{artstr_cmstr}_details.png", 
            m=["modelA", "modelB1-10", "modelB1-100"] +
             ["modelA-strongVM", "modelB2-10", "modelB2-100"] + ["LargePVS", "LargePVSA"],
             tstr=["1+2+3+4", "2+4+6+8"],
             artstr_cmstr=["MCA-R_0.2+5.0+20.0+20.0", "MCA-L_0.2+5.0+20.0+20.0",
                            "BA_3.0+40.0+40.0+10.0","ACA-A3_0.2+2.0+12.0+10.0",
                            "ACA-A2_0.2+4.0+12.0+10.0","MCA2-L_0.2+8.0+8.0+5.0",
                            "MCA2-R_0.2+8.0+8.0+5.0"])

rule runSimuation:
    conda:"environment.yml"
    input:
        volmesh=lambda wildcards: getconfig(wildcards.modelname, "mesh"),
        artmesh="mesh/networks/arteries_smooth.vtk",
        venmesh="mesh/networks/venes_smooth.vtk",
        config="configfiles/{modelname}.yml",
        csf_velocity_file=lambda wildcards: getconfig(wildcards.modelname, "csf_velocity_file"),
        arterial_velocity_file=lambda wildcards: getconfig(wildcards.modelname, "arterial_velocity_file"),
        venous_velocity_file=lambda wildcards: getconfig(wildcards.modelname, "venous_velocity_file"),
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
        #bg="mesh/standard/standard.xdmf"
    output:
        csf_dispersion_file="results/csf_flow/{csf_flow_model}/R.xdmf",
        csf_dispersion_plot="results/csf_flow/{csf_flow_model}/R.png",
    shell:
        "python scripts/compute_dispersion_field.py {input.csf_pressure_file} {output.csf_dispersion_file}"

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
        "python scripts/generate_plot.py {wildcards.modelname} {wildcards.time} {wildcards.type} --filename {output.plot}"

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
        "python scripts/create_diff_plot.py {wildcards.model1} {wildcards.model2} {wildcards.time} {wildcards.type} --filename {output.plot}"

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
        ven2="results/{model2}/{model2}_vein.xdmf",
        metrics_yaml1="results/{model1}/mean_concentrations.yml",
        metrics_yaml2="results/{model2}/mean_concentrations.yml",
    output:
        plot="plots/comparisons/{model1}_{model2}/{model1}_{model2}_{type}.png"
    params:
        cmax= lambda wildcards: cmax[wildcards.type][f"{wildcards.model1}_{wildcards.model2}"],
        diffmax= lambda wildcards: diffmax[wildcards.type][f"{wildcards.model1}_{wildcards.model2}"]
    shell:
        """
        python scripts/compare_models.py {wildcards.model1} {wildcards.model2} {wildcards.type} {params.cmax} {params.diffmax} {times} &&
        python scripts/compare_models_horizontal.py {wildcards.model1} {wildcards.model2} {wildcards.type} {params.cmax} {times} &&
        python scripts/generate_comparison_barplots.py {wildcards.model1} {wildcards.model2}
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
        "python scripts/analyze_tracer_dist.py {wildcards.modelname} {times}"

rule totalTracer:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.xdmf",
        art="results/{modelname}/{modelname}_artery.xdmf",
        ven="results/{modelname}/{modelname}_vein.xdmf",
    output:
        plot="plots/{modelname}/{modelname}_total_conc.png",
        metrics_yaml="results/{modelname}/mean_concentrations.yml",
    shell:
        "python scripts/mean_concentrations.py {wildcards.modelname} && " +
        "python scripts/plot_mean_concentrations.py {wildcards.modelname}"

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
        "data/T1_synthseg.nii.gz",
        config="configfiles/meshconfig/{meshname}.yml"
    output:
        [f"mesh/{{meshname}}/surfaces/{n}.ply" for n in ["LV", "parenchyma_incl_ventr", "skull", "V34"]]
    shell:
        """
        python scripts/extract_synthseg_surfaces.py {input.config}
        """

rule generateIdealizedSurfaces:
    conda:"environment.yml"
    output:
        [f"mesh/idealized/surfaces/{n}.ply" for n in ["LV", "parenchyma_incl_ventr", "skull", "V34"]]
    shell:
        """
        python scripts/generate_idealized_surfaces.py
        """

if config.get("meshing", True):
    rule generateMesh:
        conda:"mesh_environment.yml"
        input:
            [f"mesh/{{meshname}}/surfaces/{n}.ply" for n in ["LV", "parenchyma_incl_ventr", "skull", "V34"]],
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
        surf= lambda wildcards: "mesh/{meshname}/surfaces/parenchyma_incl_ventr.ply" if getconfig(f"meshconfig/{wildcards.meshname}", "refine") else []
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
    params:
        mpirun=lambda wildcards: "" if "idealized" in wildcards.csf_flow_model or "LowRes" in wildcards.csf_flow_model else "mpiexec -n 16"
    shell:
        """
        export OMP_NUM_THREADS={threads} && \
        {params.mpirun} python scripts/sas_flow_Hdiv.py \
        configfiles/{wildcards.csf_flow_model}.yml && \
        python3 scripts/plot_csf_flow.py results/csf_flow/{wildcards.csf_flow_model}
        """


rule computeProdPVSFlow:
    conda:"environment.yml"
    input:
        "results/csf_flow/{csf_flow_model}/flow.hdf",
    output:
        hdf='results/pvs_flow_prod/{csf_flow_model}-{network}-{radius}/pvs_flow.hdf',
        outputdir=directory("results/pvs_flow_prod/{csf_flow_model}-{network}-{radius}"),
    shell:
        """
        python scripts/pvs_flow_prod.py {wildcards.csf_flow_model} {wildcards.network} {wildcards.radius}
        python scripts/evaluate_pvs_flow.py {output.hdf}
        """


rule computePeristalticPVSFlow:
    conda:"environment.yml"
    input:
        "mesh/networks/arteries_smooth.vtk",
        "configfiles/{pvs_flow_model}.yml",
    output:
        hdf='results/pvs_flow_peristaltic/{pvs_flow_model}/pvs_flow.hdf',
        outputdir=directory("results/pvs_flow_peristaltic/{pvs_flow_model}"),
    params: 
        frequency=lambda wildcards: getconfig(wildcards.pvs_flow_model, "frequency"),
        amplitude=lambda wildcards: getconfig(wildcards.pvs_flow_model, "amplitude"),
        wavelength=lambda wildcards: getconfig(wildcards.pvs_flow_model, "wavelength"),
        beta=lambda wildcards: getconfig(wildcards.pvs_flow_model, "beta"),
    shell:
        """python scripts/peristalticflow.py --recompute \
         --frequency {params.frequency} \
         --amplitude {params.amplitude} \
         --wavelength {params.wavelength} \
         --beta {params.beta} \
         --output {output.outputdir} && \
         python scripts/evaluate_pvs_flow.py {output.hdf}
         """

rule evaluatePVSFlow:
    conda:"environment.yml"
    input:
        hdf='results/{flow_type}/{pvs_flow_model}-{netw}/pvs_flow.hdf',
    output:
        plots="plots/{flow_type}/{pvs_flow_model}-{netw}/{pvs_flow_model}-{netw}_velocity.png"
    shell:
        "python scripts/evaluate_pvs_flow.py {input.hdf}"


rule generateT1OverviewPlot:
    conda:"environment.yml"
    input:
        sas1="results/{m}/{m}_sas.xdmf",
        art1="results/{m}/{m}_artery.xdmf",
        ven1="results/{m}/{m}_vein.xdmf",    
    output:
        "plots/{m}/{m}_overview_{times}.png"
    params:
        times=lambda wildcards: wildcards.times.split("-")
    shell:
        "python scripts/overview_plot.py {wildcards.m} {params.times}"

rule detailPlot:
    conda:"environment.yml"
    input:
        sas1="results/{m}/{m}_sas.xdmf",
        art1="results/{m}/{m}_artery.xdmf",
        ven1="results/{m}/{m}_vein.xdmf",    
    output:
        "plots/{m}/{m}_{tstr}_{artstr}_{cmstr}_details.png",
    wildcard_constraints: 
        artstr=".*",
        cmstr=".*"
    params:
        times=lambda wildcards: wildcards.tstr.split("+"),
        artlabels=lambda wildcards: wildcards.artstr.split("+"),
        cmax=lambda wildcards: wildcards.cmstr.split("+")
    shell:
        "python scripts/detail_plot.py {wildcards.m} " +
        " --times  {params.times} --artlabels {params.artlabels}" +
        " --cmax {params.cmax}"


rule makeVideo:
    conda:"environment.yml"
    input:
        sas="results/{m}/{m}_sas.xdmf",
        art="results/{m}/{m}_artery.xdmf",
    output:
        "plots/{m}/{m}.mp4"
    shell:
        "python scripts/make_video.py {wildcards.m}"










    







