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
from matplotlib.markers import MarkerStyle
import yaml
from branch_marking import color_branches
from label_arteries import pointlabels
from test_map_on_global_coords_shift import map_kdtree

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

def compare_concentrations(modelname:str):
    config = read_config(f"configfiles/{modelname}.yml")
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, 3*dt*config["output_frequency"])
    sas_conc, subd_marker = get_result_fenics(modelname, "sas", times)
    art_conc, art_radii = get_result_fenics(modelname, "artery", times)
    ven_conc, ven_radii = get_result_fenics(modelname, "vein", times)
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

    #inflow = np.cumsum([k*max(0.0, t0 - t) for t in alltimes]) * dt
    set_plotting_defaults()
    sns.set_palette("magma_r", n_colors=4)
    mol2mmol = 1e3
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    tot = np.zeros_like(times)
    labels = ["CSF", "parenchyma", "PVS artery", "PVS vein"]
    compartments = [csf_tot, par_tot, art_tot, ven_tot]
    volumes = [assemble(phi*dx(CSFID)), assemble(phi*dx(PARID)),
               assemble(area_artery*dxa),assemble(area_vein*dxv)]
    for q,v, l in zip(compartments, volumes, labels):
        new_tot = tot + q*mol2mmol
        #plt.plot(times / (60*60), new_tot, ".-", label=l)
        fill = ax.fill_between(times/ (60*60), new_tot, tot, alpha=0.7, label=l)
        ax2.plot(times/ (60*60), q/v,"o-" , ms=10,
                 mfc=fill.get_facecolor(), color="grey")
        tot = new_tot
    
    #plt.plot(alltimes / (60*60), inflow, "--" , label="inflow")
    plt.xlabel("time (h)")
    ax2.set_ylabel('mean tracer concentration (mol/l)')
    ax.set_ylabel("total tracer content (mmol)")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=6, columnspacing=0.7)
    plt.tight_layout()
    filename = f"plots/{modelname}/{modelname}_total_conc.png"
    plt.savefig(filename)

    metrics = dict()
    for l,v, comp in zip(labels, volumes, compartments):
        metrics.update({f"cmean_{l}_{t}":float(c/v) for t,c in zip(times, comp)})
        metrics.update({f"ctot_{l}_{t}":float(c) for t,c in zip(times, comp)})

    with open(f'results/{modelname}/mean_concentrations.yml', 'w') as outfile:
        yaml.dump(metrics, outfile, default_flow_style=False)

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
    V = art_conc[0].function_space()
    pvs_radii = Function(DG0a)
    pvs_radii.vector().set_local(pvs_ratio_artery*art_radii.vector().get_local())
    priority = InterfaceResolution(subdomains=subd_marker,
                                       resolve_conflicts={(CSFID, PARID): 4})
    pvs_shape = Circle(radius=pvs_radii, degree=40, quad_rule='midpoint')
    proximity_dist = 2e-3
    nearby_shape = Circle(radius=project(pvs_radii + proximity_dist, DG0a), degree=40, quad_rule='midpoint')
    c_averages = [ii_project(Average(c, artery, pvs_shape, 
                                     normalize=True), V) for c in sas_conc]
    nearby_averages = [ii_project(Average(c, artery, nearby_shape, 
                                     normalize=True), V) for c in sas_conc]
    
    segments, segids ,_ = color_branches(artery)
    dxs = Measure("dx", artery, subdomain_data=segments)
    points = np.array([c for n,c in pointlabels])
    artlabels = [n for n,c in pointlabels]
    cellidx = map_kdtree(DG0a.tabulate_dof_coordinates(), points)
    artsegids = [int(segments.array()[i]) for i in cellidx]

    seglengths = [assemble(1*dxs(si)) for si in artsegids]
    conc_at_point, avg_conc_around_point, avg_conc_nearby_point= {}, {}, {}
    for n, si,l in zip(artlabels, artsegids, seglengths):
        conc_at_point[n] = [assemble(c*dxs(si))/l for c in art_conc]
        avg_conc_around_point[n] = [assemble(c*dxs(si))/l for c in c_averages]
        avg_conc_nearby_point[n] = [assemble(c*dxs(si))/l for c in nearby_averages]
    artlabelsets = [["BA",],["ICA-L", "ICA-R"], ["MCA-L", "MCA-R"],["MCA2-L", "MCA2-R"], 
                    ["ACA-A1-L", "ACA-A1-R"], ["ACA-A2", "ACA-A3"],
                    ["PCA-L", "PCA-R"], ["PER-L", "PER-R"]]
    coldict = dict(zip(artlabels, sns.color_palette("dark", n_colors=len(artlabels))))
    fig, axes = plt.subplots(2, int(np.ceil(len(artlabelsets)/2)), figsize=(14,8), constrained_layout=True)
    for ax, al_set in zip(axes.flatten(), artlabelsets):
        peaks, lags, ftas = [],[],[]
        for n in al_set:
            col = coldict[n]
            h1 = ax.plot(times / 3600, conc_at_point[n], color=col, label=n)
            h2 = ax.plot(times / 3600, avg_conc_around_point[n], color=col, ls="dashed")
            h3 = ax.plot(times / 3600, avg_conc_nearby_point[n], color=col, ls="dotted")
            ax.fill_between(times / 3600, conc_at_point[n], avg_conc_around_point[n],
                             alpha=0.2, color=col)
            ax.fill_between(times / 3600, avg_conc_around_point[n], avg_conc_nearby_point[n],
                             alpha=0.1, color=col)
            pvs_peak = times[np.argmax(conc_at_point[n])]
            avg_peak = times[np.argmax(avg_conc_around_point[n])]
            lag = avg_peak - pvs_peak
        
            fta = times[np.where(np.array(conc_at_point[n]) > 0.1)[0][0]]
            #ax.axvline(fta/3600,ls="-", ymax=0.2, color=col)
            peaks.append(pvs_peak); lags.append(lag), ftas.append(fta)
        format_secs = lambda s: "{:1.0f}:{:02.0f}".format(*divmod(abs(s)/60, 60))
        #format_secs = lambda s: f"{s/60} min"
        nl = chr(10)
        text = (f"peak: {'/ ' + (chr(10)).join(format_secs(p) for p in peaks)} h" + nl +
                f"Δt: {nl.join(('+' if l>-60 else '-') + format_secs(l) for l in lags)} h" + nl +
                f"FTA: {nl.join(format_secs(p) for p in ftas)} h")
        x_text = 1 if np.mean(peaks) < T/2 else 0.4
        ax.set_xlim((0, 12))
        ax.text(x_text + (max(peaks) > 10*60*60)*0.05, 0.9,
                text, transform=ax.transAxes, horizontalalignment='right',
                verticalalignment="top")

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, 
                  columnspacing=0.3, frameon=False)
        ax.set_xlabel("time (h)")
        ax.set_ylabel("concentration (mol/l)")
    plt.figlegend(handles =[h1[0], h2[0], h3[0]],
                  labels=["PVS", "outer PVS boundary", 
                          f"PVS proximity (R + {int(proximity_dist*1e3)} mm)"],
                  loc='lower center', #facecolor="white", 
                  edgecolor="black", frameon=False,
                  ncols=3, bbox_to_anchor=(0.5, 0.98))
    plt.savefig(f"plots/{modelname}/{modelname}_conc_at_label.png",
                bbox_inches='tight')

    timespoints = np.array([1,3,6,12,24])*3600
    timeidx = np.where(np.isin(times, timespoints))[0]
    times, c_averages, art_conc = times[timeidx], [c_averages[i] for i in timeidx], [art_conc[i] for i in timeidx]

    conc_jump = [project((cavg - cart), DG0a) for  cavg, cart in zip(c_averages, art_conc)]
    rel_conc_jump = [project(100*(cavg - cart)/(cavg  + 1e-16), DG0a) for cavg, cart in zip(c_averages, art_conc)]

    pvdjump = File(f"results/{modelname}/{modelname}_rel_jump_artery.pvd") 
    for rcj,t in zip(rel_conc_jump, times):
        rcj.rename("rel c jump", "rel c jump")
        pvdjump.write(rcj, t)

    # mmol/l equals mol/m^3
    xlabels = {"abs": "$\Delta c$ (mmol/l)", "rel":"$\Delta c_{rel}$ (%)"}
    kdecut = {"abs": 0, "rel":0}
    perclevel = {"abs": 80, "rel":95}
    for mode, jump_data in zip(["abs", "rel"], [conc_jump, rel_conc_jump]):
        cj_arr = [cj.vector() for cj in jump_data]
        jumpminmax = minmax(cj_arr, percentile=perclevel[mode])

        for cj,t in zip(jump_data, times):
            plt.figure()
            plt.hist(cj.vector()[:], range=jumpminmax, density=True, bins=100)
            filename = f"plots/{modelname}/{modelname}_{mode}_conc_jump_{t}.png"
            plt.savefig(filename)

        plt.figure()
        plt.hist(cj_arr, range=jumpminmax, density=True, bins=20)
        filename = f"plots/{modelname}/{modelname}_{mode}_conc_jump.png"
        plt.savefig(filename)

        cjdf = pd.DataFrame({f"{int(t/3600)} h":cj for t,cj in zip(times, cj_arr)})
        cjdfclipped = cjdf.copy()
        jumpminmax = (jumpminmax[0], 100)
        cjdfclipped[cjdf > jumpminmax[1]] = np.nan
        cjdfclipped[cjdf < jumpminmax[0]] = np.nan

        plt.figure()
        sns.kdeplot(cjdfclipped, bw_adjust=1, cut=kdecut[mode], fill=True, common_norm=False,
                    palette="magma", alpha=.3, linewidth=1,)
        if mode=="abs":
            plt.xlim(np.array(jumpminmax)*0.5)
        plt.xlim((-50, 100))
        filename = f"plots/{modelname}/{modelname}_{mode}_conc_jump_kde.png"
        plt.xlabel(xlabels[mode])
        plt.tight_layout()
        plt.savefig(filename)

        plt.figure()
        sns.histplot(cjdfclipped, fill=True, common_norm=False, bins=50,
                    palette="crest", alpha=.5, linewidth=1,multiple="dodge",)
        plt.xlabel(xlabels[mode])
        filename = f"plots/{modelname}/{modelname}_{mode}_conc_jump_hist.png"

        plt.tight_layout()
        plt.savefig(filename)



if __name__ == "__main__":
    typer.run(compare_concentrations)