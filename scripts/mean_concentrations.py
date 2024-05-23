from fenics import *
from xii import *
from plotting_utils import get_result_fenics
import typer
from typing import List
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config, minmax
from solver import pcws_constant, read_vtk_network, as_P0_function
import numpy as np
import seaborn as sns
from IPython import embed
import pandas as pd

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5


def compare_concentrations(modelname:str, times:List[int]):
    sas_conc, subd_marker = get_result_fenics(modelname, "sas", times)
    art_conc, art_radii = get_result_fenics(modelname, "artery", times)
    ven_conc, ven_radii = get_result_fenics(modelname, "vein", times)
    config = read_config(f"configfiles/{modelname}.yml")
    pvs_ratio_artery = config["pvs_ratio_artery"]
    pvs_ratio_vein = config["pvs_ratio_venes"]
    subd_marker.array()[np.isin(subd_marker.array(),[LVID, V34ID, CSFNOFLOWID])] = CSFID
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
    dt , T= config["dt"], config["T"]
    alltimes = np.arange(0, times[-1], dt)
    #inflow = np.cumsum([k*max(0.0, t0 - t) for t in alltimes]) * dt
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
    plt.ylabel("total tracer content (mmol)")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=6, columnspacing=0.7)
    plt.tight_layout()
    filename = f"plots/{modelname}/{modelname}_total_conc.png"
    plt.savefig(filename)

    def ii_project(f, V):
        u = TrialFunction(V)
        v = TestFunction(V)  
        dx_a = Measure('dx', domain=V.mesh())
        a  = inner(u,v)*dx_a  
        L  = inner(f, v)*dx_a
        A,b = map(ii_assemble, (a,L)) 
        A,b = map(ii_convert, (A,b)) 
        sol = Function(V) 
        solve(A, sol.vector(), b) 
        return sol

    artery = art_radii.function_space().mesh()
    DG0a = FunctionSpace(artery, "DG", 0)
    pvs_radii = Function(DG0a)
    pvs_radii.vector().set_local(pvs_ratio_artery*art_radii.vector().get_local())
    pvs_shape = Circle(radius=pvs_radii, degree=20,)
    c_averages = [ii_project(Average(c, artery, pvs_shape), DG0a) for c in sas_conc]
    conc_jump = [project(cavg - cart, DG0a) for  cavg, cart in zip(c_averages, art_conc)]
    rel_conc_jump = [project(100*(cavg - cart)/( 0.5*(cavg + cart) + 1e-16), DG0a) for cavg, cart in zip(c_averages, art_conc)]

    pvdjump = File(f"results/{modelname}/{modelname}_rel_jump_artery.pvd") 
    for rcj,t in zip(rel_conc_jump, times):
        rcj.rename("rel c jump", "rel c jump")
        pvdjump.write(rcj, t)

    xlabels = {"abs": "$\Delta c$ (mmol/l)", "rel":"$\Delta c_{rel}$ (%)"}
    kdecut = {"abs": 0, "rel":2}
    perclevel = {"abs": 80, "rel":90}
    for mode, jump_data in zip(["abs", "rel"], [conc_jump, rel_conc_jump]):
        cj_arr = [cj.vector() for cj in jump_data]
        jumpminmax = minmax(cj_arr, percentile=perclevel[mode])

        for cj,t in zip(jump_data, times):
            plt.figure()
            plt.hist(cj.vector()[:], range=jumpminmax, density=True, bins=20)
            filename = f"plots/{modelname}/{modelname}_{mode}_conc_jump_{t}.png"
            plt.savefig(filename)

        plt.figure()
        plt.hist(cj_arr, range=jumpminmax, density=True, bins=20)
        filename = f"plots/{modelname}/{modelname}_{mode}_conc_jump.png"
        plt.savefig(filename)

        cjdf = pd.DataFrame({f"{int(t/3600)} h":cj for t,cj in zip(times, cj_arr)})
        cjdfclipped = cjdf.copy()
        cjdfclipped[cjdf > jumpminmax[1]] = np.nan
        cjdfclipped[cjdf < jumpminmax[0]] = np.nan

        plt.figure()
        sns.kdeplot(cjdfclipped, bw_adjust=3, cut=kdecut[mode], fill=True, common_norm=False,
                    palette="crest", alpha=.3, linewidth=1)
        filename = f"plots/{modelname}/{modelname}_{mode}_conc_jump_kde.png"
        plt.xlabel(xlabels[mode])
        plt.tight_layout()
        plt.savefig(filename)

        plt.figure()
        sns.histplot(cjdfclipped, fill=True, common_norm=False, bins=10,
                    palette="crest", alpha=.5, linewidth=1,multiple="dodge",)
        plt.xlabel(xlabels[mode])
        filename = f"plots/{modelname}/{modelname}_{mode}_conc_jump_hist.png"

        plt.tight_layout()
        plt.savefig(filename)



if __name__ == "__main__":
    typer.run(compare_concentrations)