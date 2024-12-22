from plotting_utils import read_config, set_plotting_defaults
import typer
from typing import List
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

m2mum = 1e6

def compare_pvs_flow(flow_model:List[str]):

    v = "_".join([f.replace("/","+") for f in flow_model])
    os.makedirs(f"plots/comparisons/{v}/", exist_ok=True)

    results = dict()
    for m in flow_model:
        print(f"reading {m}")
        md = dict()
        md.update(read_config(f"plots/{m}/pvs_metrics.yml"))
        results[m] = md
    print(results)
    df = pd.DataFrame(results).transpose() * m2mum
    df["modelstr"] = df.index.astype("str")
    df["enlarged"] = "control"
    df["enlarged"][df.modelstr.str.contains("large") + df.modelstr.str.contains("3")] = "dilated"
    df["type"] = "none"
    df["type"][df.modelstr.str.contains("vasomotion") ] = "vasomotion"
    df["type"][df.modelstr.str.contains("cardiac") ] = "cardiac"
    df["type"][df.modelstr.str.contains("sas") ] = "pressure-driven"

    df["vm"] = df.modelstr.str.contains("vasomotion") 
    df["cardiac"] = df.modelstr.str.contains("cardiac") 
    df["pressure"] = df.modelstr.str.contains("sas") 

    set_plotting_defaults()
    sns.set_palette(["#2ec4b6", "#124559"])
    for var, ylabel in zip(["pvsumax", "pvsumean"],
                           ["max PVS velocity (μm/s)", "mean PVS velocity (μm/s)"]):
        plt.figure(figsize=(4.0,2.3))
        ax = sns.barplot(df, y=var, x="type", hue="enlarged", order=["pressure-driven", "cardiac", "vasomotion"])
        for p in ax.patches:
            h, w, x = p.get_height(), p.get_width(), p.get_x()
            if h>0:
                xy = (x + w / 2., h)
                text = f'{h:0.2f}'
                plt.annotate(text=text, xy=xy, ha='center', va='bottom')
        plt.ylabel(ylabel)
        plt.xlabel("")
        ax.legend_.set_title("")
        ax.legend_.set_frame_on(False)
        plt.savefig(f"plots/comparisons/{v}/{var}_comp.png", transparent=True, 
                    dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    typer.run(compare_pvs_flow)

# python3 scripts/compare_pvs_flow.py pvs_flow_peristaltic/cardiac_pvs_oscillation pvs_flow_peristaltic/cardiac_pvs_oscillation_enlarged pvs_flow_peristaltic/vasomotion-strong pvs_flow_peristaltic/vasomotion-strong-enlarged pvs_flow_prod/sas_flow-arteries pvs_flow_prod/sas_flow-arteries-3
