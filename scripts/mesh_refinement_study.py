import yaml
from plotting_utils import read_config
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

os.makedirs("plots/mesh_refinement/", exist_ok=True)

meshes = ["LowRes", "standard", "HighRes"]
csf_flow_model = {"LowRes":"sas_flow_LowRes",
                  "standard":"sas_flow",
                  "HighRes":"sas_flow_HighRes"}
transport_model = {"LowRes":"modelALowRes",
                  "standard":"modelA",
                  "HighRes":"modelAHighRes"}
results = dict()
for m in meshes:
    md = dict()
    # get mesh stats
    md.update(read_config(f"mesh/{m}/{m}.yml"))
    # get CSF flow stats
    csfm = csf_flow_model[m]
    md.update(read_config(f"results/csf_flow/{csfm}/metrics.yml"))
    # get PVS flow stats
    md.update(read_config(f"results/csf_flow/{csfm}/pvs_metrics.yml"))
    # get cardiac PVS flow stats
    cardiac_csf = read_config(f"results/csf_flow/cardiac_{csfm}/metrics.yml")
    md.update({f"cardiac_{k}":v for k,v in cardiac_csf.items()})
    # get mean concentrations
    tm = transport_model[m]
    md.update(read_config(f"results/{tm}/mean_concentrations.yml"))
    trmetrics = read_config(f"results/{tm}/{tm}_metrics.yml")
    md.update({"Rmax":trmetrics["R_max"],"Rmean":trmetrics["R_mean"],
               "Rmin":trmetrics["R_min"] })
    for dom in ["csf", "art", "ven"]:
        md.update({f"{dom}_max": max(trmetrics[f"{dom}_max"]),
                   f"{dom}_min": max(trmetrics[f"{dom}_min"]),
                   f"rel_undershot_{dom}": min(np.array(trmetrics[f"{dom}_min"]) /
                                               np.array(trmetrics[f"{dom}_max"]) ) },
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
fig.savefig("plots/mesh_refinement/table.png")


for m in meshes:
    df[f"{m}_rel"] = df[m] / df["standard"]

dfrel = df[[f"{m}_rel" for m in meshes]]

plt.figure()
#dfrel = dfrel[~dfrel.index.isin(["ven_min", "rel_undershot_ven"])]
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
plt.savefig("plots/mesh_refinement/mesh_bars.png")

# investigate undershoots
plt.figure()

for m in meshes:
    tm = transport_model[m]
    tmconfig = read_config(f"configfiles/{tm}.yml")
    trmetrics = read_config(f"results/{tm}/{tm}_metrics.yml")
    csf_max = trmetrics["csf_max"]
    csf_min = trmetrics["csf_min"]
    times = np.arange(0, tmconfig["T"], tmconfig["dt"]) / (60*60)
    plt.plot(times, csf_max, label=f"max {m}")
    plt.plot(times, csf_min, label=f"min {m}")
plt.xlabel("time (h)")
plt.ylabel("tracer concentration (mol/L)")
plt.legend()
plt.savefig(f"plots/mesh_refinement/undershoots.png")


