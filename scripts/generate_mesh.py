import wildmeshing as wm
import pyvista as pv
import numpy as np
import meshio
import json
import os
import typer
from plotting_utils import read_config
from pathlib import Path

def generate_mesh(configfile : str):

    config = read_config(configfile)
    meshname = Path(configfile).stem

    csg_tree = {"operation":"union",
                "left":f"mesh/{meshname}/surfaces/skull.ply",
                "right":f"mesh/{meshname}/surfaces/parenchyma.ply",
                #"right":{"operation":"union",
                #    "left":"mesh/surfaces/parenchyma.ply",
                #    "right":"mesh/surfaces/gm.ply",
                #        } 
                } 

    tetra = wm.Tetrahedralizer(epsilon=config["mesh_eps"], edge_length_r=config["edge_length"],
                               coarsen=True)
    tetra.load_csg_tree(json.dumps(csg_tree))
    tetra.tetrahedralize()
    point_array, cell_array, marker = tetra.get_tet_mesh()
    mesh = meshio.Mesh(
            point_array, [("tetra", cell_array)], cell_data={"label": [marker.ravel()]}
        )
    mesh = pv.from_meshio(mesh).clean()
    mesh['orig_indices'] = np.arange(mesh.n_cells, dtype=np.int32)
    sas = mesh.extract_cells(mesh["label"]==1).connectivity()
    mesh["sas_components"] = np.zeros(mesh.n_cells)
    mesh["sas_components"][sas["orig_indices"]] = sas["RegionId"] + 1
    os.makedirs(f"mesh/{meshname}/volmesh", exist_ok=True)
    pv.save_meshio(f"mesh/{meshname}/volmesh/mesh.xdmf", mesh)


if __name__ == "__main__":
    typer.run(generate_mesh)
