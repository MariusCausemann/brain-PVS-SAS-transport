import numpy as np
models = ["modelA", "modelB", "modelC", "modelE"] #, "modelD"]
times = list(np.array([1, 6, 12, 18, 24])*60*60)
conctimes =  list(np.array([0, 1,2,3, 4, 5, 6, 9, 12, 15, 18, 21, 24])*60*60)

cmax = {"detail":{"modelA_modelB":5, "modelB_modelC":5},
        "overview":{"modelA_modelB":10, "modelB_modelC":8, "modelB_modelE":10},  
        "isosurf":{"modelA_modelB":10, "modelB_modelC":8, "modelB_modelE":8,},          
        
}

diffmax = {"detail":{"modelA_modelB":1, "modelB_modelC":1},
           "overview":{"modelA_modelB":1, "modelB_modelC":5, "modelB_modelE":1},     
            "isosurf":{"modelA_modelB":1, "modelB_modelC":5, "modelB_modelE":1},               
}
types = ["overview","detail", "isosurf"]

rule all:
    input:
        "plots/comparisons/modelA_modelB/modelA_modelB_overview.png",
        "plots/comparisons/modelB_modelE/modelB_modelE_overview.png",
        "plots/comparisons/modelB_modelE/modelB_modelE_isosurf.png",
        "plots/comparisons/modelA_modelB/modelA_modelB_detail.png",
        "plots/comparisons/modelB_modelC/modelB_modelC_overview.png",
        "plots/comparisons/modelB_modelC/modelB_modelC_detail.png",
        #"plots/comparisons/modelA_modelD/modelA_modelD.png",
        expand("plots/{modelname}/{modelname}_tracer_vessel_dist.png", modelname=models),
        expand("plots/{modelname}/{modelname}_total_conc.png", modelname=models),
        expand("plots/{modelname}/{modelname}_{tp}_{t}.png", modelname=models, t=times, tp=types)


rule runSimuation:
    conda:"environment.yml"
    input:
        volmesh="mesh/T1/volmesh/mesh.xdmf",
        artmesh="mesh/networks/arteries_smooth.vtk",
        venmesh="mesh/networks/venes_smooth.vtk",
        flowfield="results/pvs_flow/uhat_prod.xdmf",
        config="configfiles/{modelname}.yml"
    output:
        sas="results/{modelname}/{modelname}_sas.pvd",
        art="results/{modelname}/{modelname}_arteries.pvd",
        ven="results/{modelname}/{modelname}_venes.pvd",
    shell:
        "python scripts/time_dependent.py {input.config}"

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
        sas="results/{modelname}/{modelname}_sas.pvd",
        art="results/{modelname}/{modelname}_arteries.pvd",
        ven="results/{modelname}/{modelname}_venes.pvd",
    output:
        plot="plots/{modelname}/{modelname}_{type}_{time,[0-9]*}.png"
    shell:
        "python scripts/generate_plot.py {wildcards.modelname} {wildcards.time} {wildcards.type} --filename {output.plot}"

rule generateDiffPlot:
    conda:"environment.yml"
    input:
        sas1="results/{model1}/{model1}_sas.pvd",
        art1="results/{model1}/{model1}_arteries.pvd",
        ven1="results/{model1}/{model1}_venes.pvd",
        sas2="results/{model2}/{model2}_sas.pvd",
        art2="results/{model2}/{model2}_arteries.pvd",
        ven2="results/{model2}/{model2}_venes.pvd",
    output:
        plot="plots/comparisons/{model1}_{model2}/{model1}_{model2}_diff_{type}_{time}.png"
    shell:
        "python scripts/create_diff_plot.py {wildcards.model1} {wildcards.model2} {wildcards.time} {wildcards.type} --filename {output.plot}"

rule getConcentrationRange:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.pvd",
        art="results/{modelname}/{modelname}_arteries.pvd",
        ven="results/{modelname}/{modelname}_venes.pvd",
    output:
        ranges="plots/{modelname}/{modelname}_ranges.yml"
    params:
        times=times
    shell:
        "python scripts/compute_ranges.py {wildcards.modelname} {params.times}"


rule compareModels:
    conda:"environment.yml"
    input:
        sas1="results/{model1}/{model1}_sas.pvd",
        art1="results/{model1}/{model1}_arteries.pvd",
        ven1="results/{model1}/{model1}_venes.pvd",
        sas2="results/{model2}/{model2}_sas.pvd",
        art2="results/{model2}/{model2}_arteries.pvd",
        ven2="results/{model2}/{model2}_venes.pvd"
    output:
        plot="plots/comparisons/{model1}_{model2}/{model1}_{model2}_{type}.png"
    params:
        cmax= lambda wildcards: cmax[wildcards.type][f"{wildcards.model1}_{wildcards.model2}"],
        diffmax= lambda wildcards: diffmax[wildcards.type][f"{wildcards.model1}_{wildcards.model2}"]
    shell:
        "python scripts/compare_models.py {wildcards.model1} {wildcards.model2} {wildcards.type} {params.cmax} {params.diffmax} {times}"

rule analyzeTracerDist:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.pvd",
        art="results/{modelname}/{modelname}_arteries.pvd",
        ven="results/{modelname}/{modelname}_venes.pvd",
    output:
        plot="plots/{modelname}/{modelname}_tracer_vessel_dist.png"
    shell:
        "python scripts/analyze_tracer_dist.py {wildcards.modelname} {times}"

rule totalTracer:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.pvd",
        art="results/{modelname}/{modelname}_arteries.pvd",
        ven="results/{modelname}/{modelname}_venes.pvd",
    output:
        plot="plots/{modelname}/{modelname}_total_conc.png"
    shell:
        "python scripts/mean_concentrations.py {wildcards.modelname} {conctimes}"

rule segmentT1:
    input:
        "data/T1.nii.gz"
    output:
        "results/freesurfer/T1_synthseg.nii.gz"
    shell:
        "scripts/synthseg.sh"

rule meshT1:
    conda:"environment.yml"
    input:
        "results/freesurfer/T1_synthseg.nii.gz"
    output:
        "mesh/T1/volmesh/mesh.xdmf",
        "mesh/T1/volmesh/mesh.h5",
    shell:
        """
        python scripts/extract_synthseg_surfaces.py &&
        python scripts/generate_synthseg_mesh.py
        """

rule computeSASFlow:
    conda:"environment.yml"
    input:
        "mesh/T1/volmesh/mesh.xdmf",
        "mesh/T1/volmesh/mesh.h5",
    output:
        "results/csf_flow/T1/csf_v.xdmf",
        "results/csf_flow/T1/csf_v.h5",
    shell:
        "python scripts/sas_flow.py T1"


rule AverageSASFlow2PVS:
    conda:"environment.yml"
    input:
        "results/csf_flow/T1/csf_v.xdmf",
        "results/csf_flow/T1/csf_v.h5",
    output:
        'results/pvs_flow/uhat_prod.xdmf',
        'results/pvs_flow/uhat_prod.h5'
    shell:
        "python scripts/uhat_prod.py"





    







