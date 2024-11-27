import wildmeshing as wm
import pyvista as pv
import numpy as np
import meshio
import json
import os
import typer
from pathlib import Path
from plotting_utils import read_config
import yaml
CSFID = 1
PARID = 2
LVID = 3
V34ID = 4
CSFNOFLOWID = 5

def get_csg_tree(folder):
    tree = {"operation":"union",
                "right":
                    {"operation":"union",
                    "left":f"{folder}/LV.ply",
                    "right":f"{folder}/V34.ply",
                    },
                "left":
                    {"operation":"union",
                    "left":f"{folder}/skull.ply",
                    "right":f"{folder}/parenchyma_incl_ventr.ply",
                    },
                } 
    print(tree)
    return tree

def generate_mesh(configfile : str):

    config = read_config(configfile)
    meshname = Path(configfile).stem

    tetra = wm.Tetrahedralizer(epsilon=config["epsilon"],
                               edge_length_r=config["edge_length_r"],
                               coarsen=True, max_threads=4, stop_quality=8,
                               max_its=30)
    csg_tree = get_csg_tree(f"mesh/{meshname}/surfaces")

    tetra.load_csg_tree(json.dumps(csg_tree))
    tetra.tetrahedralize()
    point_array, cell_array, marker = tetra.get_tet_mesh()
    mesh = meshio.Mesh(
            point_array, [("tetra", cell_array)], 
            cell_data={"label": [marker.ravel()]}
        )
    mesh = pv.from_meshio(mesh).clean()
    
    # make sure all expected labels are actually there
    assert np.isin([CSFID, PARID, LVID, V34ID,], mesh["label"]).all()

    os.makedirs(f"mesh/{meshname}/volmesh", exist_ok=True)
    pv.save_meshio(f"mesh/{meshname}/volmesh/mesh.xdmf", mesh)

if __name__ == "__main__":
    typer.run(generate_mesh)
