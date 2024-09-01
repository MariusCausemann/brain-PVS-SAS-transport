import wildmeshing as wm
import pyvista as pv
import numpy as np
import meshio
import json
import os
import dolfin as df
import networkx as nx
import xii
import typer
from pathlib import Path
from plotting_utils import read_config
import yaml

CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

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

csg_tree = {"operation":"union",
            "right":
                {"operation":"union",
                 "left":"mesh/T1/surfaces/LV.ply",
                 "right":"mesh/T1/surfaces/V34.ply",
                },
            "left":
                {"operation":"union",
                 "left":"mesh/T1/surfaces/skull.ply",
                 "right":"mesh/T1/surfaces/parenchyma_incl_ventr.ply",
                },
            } 

def generate_mesh(configfile : str):

    config = read_config(configfile)
    meshname = Path(configfile).stem

    tetra = wm.Tetrahedralizer(epsilon=config["epsilon"],
                               edge_length_r=config["edge_length_r"],
                               coarsen=False)
    
    tetra.load_csg_tree(json.dumps(csg_tree))
    tetra.tetrahedralize()
    point_array, cell_array, marker = tetra.get_tet_mesh()
    mesh = meshio.Mesh(
            point_array, [("tetra", cell_array)], cell_data={"label": [marker.ravel()]}
        )
    mesh = pv.from_meshio(mesh).clean()

    # compute distance to interface for later refinement
    parenchyma_surf = pv.read("mesh/T1/surfaces/parenchyma.ply")
    dist = mesh.compute_implicit_distance(parenchyma_surf).ptc()
    mesh["par_dist"] = dist["implicit_distance"]
    os.makedirs(f"mesh/{meshname}/volmesh", exist_ok=True)
    pv.save_meshio(f"mesh/{meshname}/volmesh/mesh.xdmf", mesh)

    # read in again and use fenics to exclude problematic parts 
    # of the CSF space from the Stokes computation (disconnected domains
    # and overconstrained cells) 

    sas = df.Mesh()
    with df.XDMFFile(f'mesh/{meshname}/volmesh/mesh.xdmf') as f:
        f.read(sas)
        gdim = sas.geometric_dimension()
        label = df.MeshFunction('size_t', sas, gdim, 0)
        par_dist = df.MeshFunction('size_t', sas, gdim, 0)
        f.read(label, 'label')
        f.read(par_dist, 'par_dist')

    def refine_region(sas, label, labelids):
        mf = df.MeshFunction("bool", sas, 3, 0)
        mf.array()[:] = np.isin(label.array()[:], labelids)
        sas = df.refine(sas, mf)
        label = df.adapt(label, sas)
        return sas, label

    def refine_sphere(sas, coords, radius, label, criterion=None):
        mf = df.MeshFunction("bool", sas, 3, 0)
        rm = df.CompiledSubDomain("(x[0] - c0)*(x[0] - c0) + (x[1] - c1)*(x[1] - c1) + (x[2] - c2)*(x[2] - c2) < r*r",
                                r = 0.01, c0=coords[0],  c1=coords[1],  c2=coords[2])
        rm.mark(mf, True)
        if criterion is not None:
            mf.array()[criterion==False] = False
        sas = df.refine(sas, mf)
        label = df.adapt(label, sas)
        return sas, label

    coords = (0.0847, 0.0833, 0.001)
    crit = abs(par_dist.array()[:]) < 0.002
    sas, label = refine_sphere(sas, coords, 0.01, label, criterion=crit)
    sas, label = refine_sphere(sas, coords, 0.01, label)
    sas, label = refine_region(sas, label, labelids=[V34ID])


    coords = (0.084, 0.11, 0.052)
    sas, label = refine_sphere(sas, coords, 0.02, label)

    sas_outer = xii.EmbeddedMesh(label, [CSFID,LVID,V34ID]) 
    cell_f = color_connected_components(sas_outer) 

    colors = df.MeshFunction('size_t', sas, gdim, 0)

    cellmap = sas_outer.parent_entity_map[sas.id()][3]

    for i,m in enumerate(cell_f.array()):
        colors.array()[cellmap[i]] = m

    label.array()[colors.array()[:] >=2]  = CSFNOFLOWID

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

    label.rename("label", "label")

    with df.XDMFFile(f'mesh/{meshname}/{meshname}.xdmf') as f:
        f.write(label)
    mesh = label.mesh()
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
    typer.run(generate_mesh)
