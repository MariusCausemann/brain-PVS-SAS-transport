import wildmeshing as wm
import pyvista as pv
import numpy as np
import meshio
import json


csg_tree = {"operation":"union",
            "left":"mesh/surfaces/skull.ply",
            "right":"mesh/surfaces/gm.ply"} 

tetra = wm.Tetrahedralizer(epsilon=0.002, edge_length_r=0.05,
                           coarsen=False)
tetra.load_csg_tree(json.dumps(csg_tree))
tetra.tetrahedralize()
point_array, cell_array, marker = tetra.get_tet_mesh()
mesh = meshio.Mesh(
        point_array, [("tetra", cell_array)], cell_data={"label": [marker.ravel()]}
    )
mesh = pv.from_meshio(mesh).clean()
pv.save_meshio("mesh/volmesh/mesh.xdmf", mesh)