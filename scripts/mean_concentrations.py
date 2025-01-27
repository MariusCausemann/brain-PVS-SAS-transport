from fenics import *
from xii import *
from plotting_utils import get_result_fenics
import typer
from typing import List
import matplotlib.pyplot as plt
from plotting_utils import set_plotting_defaults, read_config, minmax
from solver import pcws_constant
import numpy as np
import seaborn as sns
from IPython import embed
import pandas as pd
import yaml
from branch_marking import color_branches
from label_arteries import pointlabels
from test_map_on_global_coords_shift import map_kdtree
from peristalticflow import mesh_to_weighted_graph, nx
from subdomain_ids import CSFID, PARID, LVID, V34ID, CSFNOFLOWID

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

def compare_concentrations(modelname:str):
    config = read_config(f"configfiles/{modelname}.yml")
    dt , T= config["dt"], config["T"]
    times = np.arange(0, T + dt, 2*dt*config["output_frequency"])
    assert 3600 in times
    sas_conc, subd_marker = get_result_fenics(modelname, "sas", times)
    art_conc, art_radii, art_roots = get_result_fenics(modelname, "artery", times, getroots=True)
    ven_conc, ven_radii, ven_roots = get_result_fenics(modelname, "vein", times, getroots=True)
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

    #inflow = np.cumsum([k*max(0.0, t0 - t) for t in alltimes]) * dt
    set_plotting_defaults()
    fig, ax = plt.subplots()
    labels = ["CSF", "parenchyma", "PVS artery", "PVS vein"]
    compartments = [csf_tot, par_tot, art_tot, ven_tot]
    volumes = [assemble(phi*dx(CSFID)), assemble(phi*dx(PARID)),
               assemble(area_artery*dxa),assemble(area_vein*dxv)]

    metrics = dict()
    for l,v, comp in zip(labels, volumes, compartments):
        metrics.update({f"cmean_{l}_{t}":float(c/v) for t,c in zip(times, comp)})
        metrics.update({f"ctot_{l}_{t}":float(c) for t,c in zip(times, comp)})


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
    csf_prio = InterfaceResolution(subdomains=subd_marker,
                                       resolve_conflicts={(CSFID, PARID): 1e6})
    c_averages = [ii_project(Average(c, artery, pvs_shape, resolve_interfaces=priority,
                                     normalize=True), V) for c in sas_conc]
    nearby_averages = [ii_project(Average(c, artery, nearby_shape, resolve_interfaces=priority,
                                     normalize=True), V) for c in sas_conc]
    c_averages_csf = [ii_project(Average(c, artery, pvs_shape, resolve_interfaces=csf_prio,
                                     normalize=True), V) for c in sas_conc]
    nearby_averages_csf = [ii_project(Average(c, artery, nearby_shape, resolve_interfaces=csf_prio,
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
        avg_conc_around_point[n] = [assemble(c*dxs(si))/l for c in c_averages_csf]
        avg_conc_nearby_point[n] = [assemble(c*dxs(si))/l for c in nearby_averages_csf]
    metrics["avg_conc_around_point"] = avg_conc_around_point
    metrics["avg_conc_nearby_point"] = avg_conc_nearby_point
    metrics["conc_at_point"] = conc_at_point
    metrics["c_averages"] = [interpolate(cavg, DG0a).vector()[:] for cavg in c_averages]
    metrics["c_pvs"] = [interpolate(cpvs, DG0a).vector()[:] for cpvs in art_conc]

    pvsftadict, avgftadict = dict(),dict() 
    for n in artlabels:
        for c, d in zip([conc_at_point,avg_conc_around_point], [pvsftadict, avgftadict] ):
            if (np.array(c[n]) > 0.1).sum() == 0: fta = None
            else:
                fta = times[np.where(np.array(c[n]) > 0.1)[0][0]]
            d[n] = fta

    pvspeaktimedict = {n:times[np.argmax(conc_at_point[n])]for n in artlabels }
    avgpeaktimedict = {n:times[np.argmax(avg_conc_around_point[n])] for n in artlabels}
    lagdict = {n: avgpeaktimedict[n] - pvspeaktimedict[n] for n in artlabels}

    for n in artlabels:
        metrics.update({f"{n}_fta":pvsftadict[n], f"{n}_avg_fta":avgftadict[n], 
                        f"{n}_lag":lagdict[n], 
                        f"{n}_pvs_peak_time":pvspeaktimedict[n],
                        f"{n}_avg_peak_time":avgpeaktimedict[n]})
        
    # compute values over distance
    G = mesh_to_weighted_graph(artery, art_radii.vector()[:])
    path_lengths = []
    for rootnode in np.where(art_roots.array()==2)[0]:
        pl = nx.single_source_dijkstra_path_length(G, rootnode, weight="length")
        path_lengths.append(nx.utils.dict_to_numpy_array(pl))
    idx = map_kdtree(artery.coordinates(), art_conc[0].function_space().tabulate_dof_coordinates())
    dist_to_root = Function(art_conc[0].function_space())
    dist_to_root.vector()[:] = np.array(path_lengths).min(axis=0)[idx]
    dist_to_root = interpolate(dist_to_root, DG0a)
    pointidx = map_kdtree(dist_to_root.function_space().tabulate_dof_coordinates(), points)
    labeldist = {n:dist_to_root.vector()[pi] for n,pi in zip(artlabels, pointidx)}
    for ln, dist in labeldist.items():
        metrics[f"{ln}_root_dist"] = dist
    cellvol = project(CellDiameter(artery)*area_artery, DG0a).vector()[:]
    metrics["cellvol"] = cellvol
    metrics["dist_to_root"] = dist_to_root.vector()[:]

    with open(f'results/{modelname}/mean_concentrations.yml', 'w') as outfile:
            yaml.dump(metrics, outfile, default_flow_style=False)


if __name__ == "__main__":
    typer.run(compare_concentrations)