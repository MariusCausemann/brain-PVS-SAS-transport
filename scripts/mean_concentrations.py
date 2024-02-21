from fenics import *
from plotting_utils import get_result_fenics
import typer
from typing import List
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config
from solver import pcws_constant
import numpy as np
import seaborn as sns
from IPython import embed

CSFID, PARID = 1,2

def compare_concentrations(modelname:str, times:List[int]):
    sas_conc, subd_marker = get_result_fenics(modelname, "sas", times)
    art_conc, art_radii = get_result_fenics(modelname, "artery", times)
    ven_conc, ven_radii = get_result_fenics(modelname, "vein", times)
    config = read_config(f"configfiles/{modelname}.yml")
    pvs_ratio_artery = config["pvs_ratio_artery"]
    pvs_ratio_vein = config["pvs_ratio_venes"]
    dx = Measure("dx", sas_conc[0].function_space().mesh(), subdomain_data=subd_marker)    
    dxa = Measure("dx", art_conc[0].function_space().mesh())
    dxv = Measure("dx", ven_conc[0].function_space().mesh())

    area_artery  = np.pi*((art_radii*pvs_ratio_artery)**2 - (art_radii)**2)
    area_vein  = np.pi*((ven_radii*pvs_ratio_vein)**2 - (ven_radii)**2)

    ecs_share = config["ecs_share"]

    phi_dict = {CSFID: Constant(1),  PARID: Constant(ecs_share)}
    phi = pcws_constant(subd_marker, phi_dict)

    par_tot = np.array([assemble(phi*c*dx(PARID)) for c in sas_conc])
    csf_tot = np.array([assemble(phi*c*dx(CSFID)) for c in sas_conc])
    art_tot = np.array([assemble(area_artery*c*dxa) for c in art_conc])
    ven_tot = np.array([assemble(area_vein*c*dxv) for c in ven_conc])
    tot = par_tot + csf_tot + art_tot + ven_tot
    times = np.array(times)
    k, t0 = 1e-9, 7200, 
    dt , T= config["dt"], config["T"]
    alltimes = np.arange(0, times[-1], dt)
    inflow = np.cumsum([k*max(0.0, t0 - t) for t in alltimes]) * dt
    set_plotting_defaults()
    sns.set_palette("dark")
    plt.figure()
    plt.plot(times / (60*60), par_tot, ".-" , label="par")
    plt.plot(times / (60*60), csf_tot, ".-" , label="csf")
    plt.plot(times / (60*60), art_tot, ".-" , label="art")
    plt.plot(times / (60*60), ven_tot, ".-" , label="ven")
    plt.plot(times / (60*60), tot, ".-" , label="total")
    #plt.plot(alltimes / (60*60), inflow, "--" , label="inflow")
    plt.legend()
    plt.xlabel("time (h)")
    plt.ylabel("total tracer")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=6, columnspacing=0.7)
    plt.tight_layout()
    filename = f"plots/{modelname}/{modelname}_total_conc.png"
    plt.savefig(filename)
    #from IPython import embed; embed()

if __name__ == "__main__":
    typer.run(compare_concentrations)