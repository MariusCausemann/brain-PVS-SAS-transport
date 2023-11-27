models = ["modelA", "modelB"]
times = [60*60*1, 60*60*6, 60*60*12, 60*60*18, 60*60*24] # secs

rule all:
    input:
        "plots/comparisons/modelA_modelB/modelA_modelB.png"

rule runSimuation:
    conda:"environment.yml"
    input:
        volmesh="mesh/volmesh/mesh.xdmf",
        artmesh="mesh/networks/arteries_smooth.vtk",
        venmesh="mesh/networks/venes_smooth.vtk",
        config="configfiles/{modelname}.yml"
    output:
        sas="results/{modelname}/{modelname}_sas.pvd",
        art="results/{modelname}/{modelname}_arteries.pvd",
        ven="results/{modelname}/{modelname}_venes.pvd",
    shell:
        "python scripts/time_dependent.py {input.config}"

rule generatePlot:
    conda:"environment.yml"
    input:
        sas="results/{modelname}/{modelname}_sas.pvd",
        art="results/{modelname}/{modelname}_arteries.pvd",
        ven="results/{modelname}/{modelname}_venes.pvd",
    output:
        plot="plots/{modelname}/{modelname}_{time}.png"
    shell:
        "python scripts/generate_plot.py {wildcards.modelname} {wildcards.time}"

rule compareModels:
    conda:"environment.yml"
    input:
        plots1=expand("plots/{{model1}}/{{model1}}_{time}.png", time=times),
        plots2=expand("plots/{{model2}}/{{model2}}_{time}.png", time=times),
    output:
        plot="plots/comparisons/{model1}_{model2}/{model1}_{model2}.png"
    shell:
        "python scripts/compare_models.py {wildcards.model1} {wildcards.model2} {times}"

