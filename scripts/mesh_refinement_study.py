import yaml
from plotting_utils import read_config
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
from plotting_utils import time_str, set_plotting_defaults


csf_flow_model = {"LowRes":"sas_flow_LowRes",
                "standard":"sas_flow",
                "HighRes":"sas_flow_HighRes",
                "large dt":"sas_flow",
                "small dt":"sas_flow"}
transport_model = {"LowRes":"modelALowRes",
                "standard":"modelA",
                "HighRes":"modelAHighRes",
                "large dt":"modelAlargedt",
                "small dt":"modelAsmalldt"}

mesh_refinement_models = ["LowRes", "standard", "HighRes"]
time_refinement_models = ["large dt", "standard", "small dt"]
variants = ["mesh_refinement", "time_refinement"]

for v, models in zip(variants, [mesh_refinement_models, time_refinement_models]):
    os.makedirs(f"plots/{v}/", exist_ok=True)


    results = dict()
    for m in models:
        md = dict()
        tm = transport_model[m]
        tmconfig = read_config(f"configfiles/{tm}.yml")
        actual_mesh = Path(tmconfig["mesh"])
        print(actual_mesh)
        # get mesh stats
        md.update(read_config(actual_mesh.with_suffix(".yml")))
        # get CSF flow stats
        csfm = csf_flow_model[m]
        md.update(read_config(f"results/csf_flow/{csfm}/metrics.yml"))
        # get PVS flow stats
        md.update(read_config(f"results/csf_flow/{csfm}/pvs_metrics.yml"))
        # get cardiac PVS flow stats
        cardiac_csf = read_config(f"results/csf_flow/cardiac_{csfm}/metrics.yml")
        md.update({f"cardiac_{k}":v for k,v in cardiac_csf.items()})
        # get mean concentrations
        md.update(read_config(f"results/{tm}/mean_concentrations.yml"))
        trmetrics = read_config(f"results/{tm}/{tm}_metrics.yml")
        md.update({"Rmax":trmetrics["R_max"],"Rmean":trmetrics["R_mean"],
                "Rmin":trmetrics["R_min"] })
        for dom in ["csf", "art", "ven"]:
            md.update({f"{dom}_max": max(trmetrics[f"{dom}_max"]),
                    f"{dom}_min": min(trmetrics[f"{dom}_min"]),
                    f"rel_undershot_{dom}": abs(min(np.array(trmetrics[f"{dom}_min"]) /
                                                np.array(trmetrics[f"{dom}_max"]) )) },
                    )
        results[m] = md

    df = pd.DataFrame(results)
    dfstr = df.applymap('{:,.3g}'.format)
    fig, ax = plt.subplots(figsize=(6,22))
    ax.axis('tight')
    ax.axis('off')
    the_table = ax.table(cellText=dfstr.values,
                            rowLabels=dfstr.index,
                            colLabels=dfstr.columns,
                            loc="center")
    fig.tight_layout()
    fig.savefig(f"plots/{v}/table.png")


    for m in models:
        df[f"{m}_rel"] = df[m] / df["standard"]

    dfrel = df[[f"{m}_rel" for m in models]]

    plt.figure()
    dfrel = dfrel[~dfrel.index.isin(["rel_undershot_csf","rel_undershot_art", "rel_undershot_ven"])]
    ax = dfrel.plot(kind="bar", figsize=(40,4), ylim=(0, 5))
    pps = ax.containers[1]
    for i, p, val in zip(range(len(pps)), pps, df["standard"].values):
        height = p.get_height()
        ax.annotate('{0:6.2g}'.format(val),
            xy=(p.get_x() + p.get_width() / 2, height),
            xytext=(0, 3 + 2*(i%2)), # 3 points vertical offset
            textcoords="offset points",
            fontsize=6,
            ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(f"plots/{v}/mesh_bars.png")

    # investigate undershoots
    set_plotting_defaults()
    linestyles = ["dotted", "dashed", "solid"]
    for dom in ["csf","art", "ven"]:
        plt.figure()

        for m,ls  in zip(models, linestyles):
            tm = transport_model[m]
            tmconfig = read_config(f"configfiles/{tm}.yml")
            trmetrics = read_config(f"results/{tm}/{tm}_metrics.yml")
            dmax = trmetrics[f"{dom}_max"]
            dmin = trmetrics[f"{dom}_min"]
            times = np.arange(0, tmconfig["T"], tmconfig["dt"]) / (60*60)
            plt.plot(times, dmax, label=f"max {m}", ls=ls, color="maroon")
            plt.plot(times, dmin, label=f"min {m}", ls=ls, color="teal")
        plt.xlabel("time (h)")
        plt.ylabel("tracer concentration (mol/L)")
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3,
                  columnspacing=0.7)
        plt.savefig(f"plots/{v}/{dom}_undershoots.png")


