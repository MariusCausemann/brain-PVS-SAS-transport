import yaml
from plotting_utils import read_config
import numpy as np

meshes = ["LowRes", ] #, "standard", "HighRes"]
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
