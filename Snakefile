models = ["modelA", "modelB"]#, "modelC", "modelD"]
times = [60*60*1, 60*60*6, 60*60*12, 60*60*18, 60*60*24] # secs

rule all:
    input:
        "plots/comparisons/modelA_modelB/modelA_modelB.png",
        "plots/comparisons/modelA_modelC/modelA_modelC.png",
        #"plots/comparisons/modelA_modelD/modelA_modelD.png",
        expand("plots/{modelname}/{modelname}_tracer_vessel_dist.png", modelname=models)

rule runSimuation:
    conda:"environment.yml"
    input:
        volmesh="mesh/volmesh/mesh.xdmf",
        artmesh="mesh/networks/arteries_smooth.vtk",
        venmesh="mesh/networks/venes_smooth.vtk",
        flowfield="results/pvs_flow/pvs_flow.xdmf",
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
        plot="plots/{modelname}/{modelname}_{time,[0-9]*}.png"
    shell:
        "python scripts/generate_plot.py {wildcards.modelname} {wildcards.time}"

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
        plot="plots/comparisons/{model1}_{model2}/{model1}_{model2}_diff_{time}.png"
    shell:
        "python scripts/create_diff_plot.py {wildcards.model1} {wildcards.model2} {wildcards.time}"

rule compareModels:
    conda:"environment.yml"
    input:
        plots1=expand("plots/{{model1}}/{{model1}}_{time}.png", time=times),
        plots2=expand("plots/{{model2}}/{{model2}}_{time}.png", time=times),
        plotsdiff=expand("plots/comparisons/{{model1}}_{{model2}}/{{model1}}_{{model2}}_diff_{time}.png", time=times),

    output:
        plot="plots/comparisons/{model1}_{model2}/{model1}_{model2}.png"
    shell:
        "python scripts/compare_models.py {wildcards.model1} {wildcards.model2} {times}"

rule analyzeTracerDist:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.pvd",
        art="results/{modelname}/{modelname}_arteries.pvd",
        ven="results/{modelname}/{modelname}_venes.pvd",
    output:
        plot="plots/{modelname}/{modelname}_tracer_vessel_dist.png"
    shell:
        "python scripts/analyze_tracer_dist.py {wildcards.modelname}"
