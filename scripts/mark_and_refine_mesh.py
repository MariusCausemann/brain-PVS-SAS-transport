import pyvista as pv
import numpy as np
import dolfin as df
import networkx as nx
import xii
import typer
from pathlib import Path
from plotting_utils import read_config
import yaml
from solver import mark_internal_interface, mark_external_boundary
from test_map_on_global_coords_shift import map_kdtree

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

# interface/boundary ids
UPPER_SKULL_ID = 1
LOWER_SKULL_ID = 2
LV_INTERF_ID = 3
PIA_ID = 4
SPINAL_CORD_ID = 5
SPINAL_OUTLET_ID = 6
CSF_NO_FLOW_CSF_ID = 7


def color_connected_components(mesh):
    '''
    Cell function with colors corresponding to connected components of mesh graph. 
    Largest components first.
    '''
    print(f'Coloring components in mesh with {mesh.num_cells()} cells')
    tdim = mesh.topology().dim()
    mesh.init(tdim-1, tdim)  #  Facet to cell

    g = nx.Graph()
    g.add_edges_from(tuple(facet.entities(tdim))
                     for facet in df.facets(mesh) if len(facet.entities(tdim)) > 1)

    ccs = sorted(list(nx.algorithms.connected_components(g)))

    
    components = df.MeshFunction('size_t', mesh, tdim, 0)
    values = components.array()
    for (color, component) in enumerate(ccs, 1):
        print(f'Component {color} with {len(component)} cells')
        values[list(component)] = color
    return components

def refine_region(sas, label, labelids):
    mf = df.MeshFunction("bool", sas, 3, 0)
    mf.array()[:] = np.isin(label.array()[:], labelids)
    sas = df.refine(sas, mf)
    label = df.adapt(label, sas)
    return sas, label

def refine_sphere(sas, coords, radius, label, criterion=None, min_size=None):
    mf = df.MeshFunction("bool", sas, 3, 0)
    rm = df.CompiledSubDomain("(x[0] - c0)*(x[0] - c0) +" +
        "(x[1] - c1)*(x[1] - c1) + (x[2] - c2)*(x[2] - c2) < r*r",
        r = radius, c0=coords[0],  c1=coords[1],  c2=coords[2])
    rm.mark(mf, True)
    if criterion is not None:
        mf.array()[criterion==False] = False
    if min_size is not None:
        cd = df.CellDiameter(sas)
        cd = df.project(cd, df.FunctionSpace(sas, "DG", 0), solver_type="cg", 
                                                        preconditioner_type="jacobi")
        mf.array()[cd.vector()[:] <= min_size] = False

    sas = df.refine(sas, mf)
    label = df.adapt(label, sas)
    return sas, label

def get_surface_dist(mesh, surface):
    DG = df.FunctionSpace(mesh, "DG", 0)
    center_cloud = pv.PolyData(DG.tabulate_dof_coordinates())
    dist = center_cloud.compute_implicit_distance(surface)
    return dist["implicit_distance"]

def mark_inlet(mesh, bm, label, vec, threshold, zmin=0):
    for f in df.facets(mesh):
        if f.exterior()==False:continue
        if f.midpoint()[2] > zmin: continue
        if np.dot(f.normal()[:], vec) > threshold:
            bm[f] = label

def mark_boundaries(label):
    mesh = label.mesh()
    bm = df.MeshFunction("size_t", mesh, 2, 0)

    # set upper and lower skull
    upper_skull_z = 0.00
    lower_skull = df.CompiledSubDomain("on_boundary")
    lower_skull.mark(bm, LOWER_SKULL_ID)
    upper_skull = df.CompiledSubDomain("on_boundary && (x[2] - 0.8*x[1]) > m",
                                     m=upper_skull_z)
    upper_skull.mark(bm, UPPER_SKULL_ID)

    # mark LV surface
    mark_internal_interface(mesh, label, bm, LV_INTERF_ID,
                            doms=[LVID, PARID])
    #mark Pia surface
    for ID in CSFID, CSFNOFLOWID, V34ID:
        mark_internal_interface(mesh, label, bm, PIA_ID,
                                doms=[ID, PARID])

    mark_internal_interface(mesh, label, bm, CSF_NO_FLOW_CSF_ID,
                            doms=[CSFID, CSFNOFLOWID])
    # mark spinal outlet
    mark_inlet(mesh, bm, SPINAL_OUTLET_ID, (0 ,0,-1), 0.8,
               zmin=mesh.coordinates()[:,2].min() + 0.5e-3)

    # mark spinal cord
    mark_external_boundary(mesh, label, bm, SPINAL_CORD_ID,
                        doms=[PARID])
    return bm

def mark_and_refine(configfile : str):

    config = read_config(configfile)
    meshname = Path(configfile).stem

    # compute distance to interface for later refinement
    parenchyma_surf = pv.read(f"mesh/{meshname}/surfaces/parenchyma_incl_ventr.ply")

    # read in again and use fenics to exclude problematic parts 
    # of the CSF space from the Stokes computation (disconnected domains)

    sas = df.Mesh()
    with df.XDMFFile(f'mesh/{meshname}/volmesh/mesh.xdmf') as f:
        f.read(sas)
        gdim = sas.geometric_dimension()
        label = df.MeshFunction('size_t', sas, gdim, 0)
        f.read(label, 'label')

    if config.get("refine", True):
        # refine V3 and V4
        sas, label = refine_region(sas, label, labelids=[V34ID])

        # refine AQ
        #AQ_coords = (0.0837, 0.08, 0.065) 
        #sas, label = refine_sphere(sas, AQ_coords, 0.007, label)

        # refine Par boundary
        print("======================")
        print("refining pia boundary...")
        bottom_refine_coords = (0.0847, 0.11, 0.001)
        for i in range(2):
            print("======================")
            crit = abs(get_surface_dist(sas, parenchyma_surf)) < 0.003
            print(sas.num_cells())
            sas, label = refine_sphere(sas, bottom_refine_coords, 0.08, label, 
                                        criterion=crit, min_size=config["edge_length_r"] / 3)
            print(sas.num_cells())

        # refine inlet area
        print("======================")
        print("refining inlet area...")
        inlet_coords = (0.0847, 0.0833, 0.001)
        for i in range(2):
            print("======================")
            print(sas.num_cells())
            crit = abs(get_surface_dist(sas, parenchyma_surf)) < 0.002
            sas, label = refine_sphere(sas, inlet_coords, 0.01, label, criterion=crit,
                                        min_size=config["edge_length_r"] / 8)
            print(sas.num_cells())

        # refine AQ area
        print("======================")
        print("refining AQ area...")
        aq_coords = (0.0815, 0.082, 0.0665)
        for i in range(1):
            print("======================")
            print(sas.num_cells())
            sas, label = refine_sphere(sas, aq_coords, 0.005, label)
            print(sas.num_cells())

    # remove unconnected components of the CSF space
    sas_outer = xii.EmbeddedMesh(label, [CSFID,LVID,V34ID])
    cell_f = color_connected_components(sas_outer) 
    colors = df.MeshFunction('size_t', sas, gdim, 0)
    cellmap = sas_outer.parent_entity_map[sas.id()][3]

    for i,m in enumerate(cell_f.array()):
        colors.array()[cellmap[i]] = m

    cols,  color_counts = np.unique(colors.array(), return_counts=True)

    label.array()[np.isin(colors.array()[:],
                          cols[np.argsort(color_counts)[:-2]])]  = CSFNOFLOWID
    label.rename("label", "label")

    def count_external_facets(mesh):
        ext_facet_count = df.MeshFunction('size_t', mesh, gdim, 0)
        tdim = mesh.topology().dim()
        mesh.init(tdim-1, tdim)  #  Facet to cell
        for i,c in enumerate(df.cells(mesh)):
            n = 0
            for f in df.facets(c):
                n+= f.exterior()
            ext_facet_count.array()[i] = n
        return ext_facet_count

    for i in range(10):
        sas_outer = xii.EmbeddedMesh(label, [CSFID,LVID,V34ID]) 
        fct_cnt = count_external_facets(sas_outer)
        marker = df.MeshFunction('bool', sas_outer, gdim, 0)
        # cells with 3 or more external facets are overconstrained
        marker.array()[:] = fct_cnt.array()[:] >= 3 
        print(f"found {marker.array().sum()} overconstrained cells")
        if marker.array().sum() == 0:
            print("no overconstrained cells, good to go!");break
        else:
            ext_marker = df.MeshFunction('bool', sas, gdim, 0)
            cellmap = sas_outer.parent_entity_map[sas.id()][3]
            for i,m in enumerate(marker.array()):
                ext_marker.array()[cellmap[i]] = m
            label.array()[ext_marker.array()] = CSFNOFLOWID 

    bm = mark_boundaries(label)

    with df.XDMFFile(f'mesh/{meshname}/{meshname}.xdmf') as f:
        f.write(label)
    with df.XDMFFile(f'mesh/{meshname}/{meshname}_facets.xdmf') as f:
        f.write(bm)

    sas_outer = xii.EmbeddedMesh(label, [CSFID,LVID,V34ID])
    bm_outer = df.MeshFunction("size_t", sas_outer, 2, 0)
    label_outer = df.MeshFunction("size_t", sas_outer, 3, 0)
    bm_ids = [UPPER_SKULL_ID,LOWER_SKULL_ID,LV_INTERF_ID,PIA_ID,
              SPINAL_CORD_ID, SPINAL_OUTLET_ID, CSF_NO_FLOW_CSF_ID]
    sas_outer.translate_markers(bm,bm_ids, marker_f=bm_outer)

    idx = map_kdtree(df.FunctionSpace(label.mesh(), "DG", 0).tabulate_dof_coordinates(),
                     df.FunctionSpace(sas_outer, "DG", 0).tabulate_dof_coordinates())
    label_outer.array()[:] = label.array()[idx]
    with df.XDMFFile(f'mesh/{meshname}/{meshname}_outer.xdmf') as f:
        f.write(label_outer)
    with df.XDMFFile(f'mesh/{meshname}/{meshname}_outer_facets.xdmf') as f:
        f.write(bm_outer)

    mesh = label.mesh()

    ds_outer = df.Measure("ds", sas_outer, subdomain_data=bm_outer)
    ds = df.Measure("ds", mesh, subdomain_data=bm)
    for i in np.unique(bm.array()[:]):
        if i == 0:
            assert df.assemble(1*ds_outer(int(i))) == 0
            assert df.assemble(1*ds(int(i))) == 0
        else:
            assert df.assemble(1*ds(int(i))) + df.assemble(1*ds_outer(int(i))) > 0

    meshinfo = {"hmax":mesh.hmax(),"hmin":mesh.hmin(),
                "ncells":mesh.num_cells(),"npoints":mesh.num_vertices()}

    labelids = {"CSFID":CSFID, "PARID":PARID, "LVID":LVID,
                "V34ID":V34ID,"CSFNOFLOWID":CSFNOFLOWID}

    dx = df.Measure("dx", mesh, subdomain_data=label)
    for k,v in labelids.items():
        meshinfo[f"vol_{k}"] = df.assemble(df.Constant(1)*dx(v))

    with open(f'mesh/{meshname}/{meshname}.yml', 'w') as outfile:
        yaml.dump(meshinfo, outfile, default_flow_style=False)

if __name__ == "__main__":
    typer.run(mark_and_refine)
