
from plotting_utils import get_result, read_config
import numpy as np
from extract_vessels import get_tubes
import typer
import pyvista as pv
import meshio
from meshio._vtk_common import vtk_to_meshio_type


def get_cells(mesh):
    if not isinstance(mesh, pv.UnstructuredGrid):
        mesh = mesh.cast_to_unstructured_grid()

    # Copy useful arrays to avoid repeated calls to properties
    vtk_offset = mesh.offset
    vtk_cells = mesh.cells
    vtk_cell_type = mesh.celltypes

    # Check that meshio supports all cell types in input mesh
    pixel_voxel = {8, 11}  # Handle pixels and voxels
    for cell_type in np.unique(vtk_cell_type):
        if cell_type not in vtk_to_meshio_type.keys() and cell_type not in pixel_voxel:
            raise TypeError(f'meshio does not support VTK type {cell_type}.')

    # Get cells
    cells = []  # type: ignore[var-annotated]
    c = 0
    for i, (offset, cell_type) in enumerate(zip(vtk_offset, vtk_cell_type)):
        if cell_type == 42:
            cell_ = mesh.get_cell(i)
            cell = [face.point_ids for face in cell_.faces]
            cell_type = f'polyhedron{cell_.n_points}'

        else:
            numnodes = vtk_cells[offset + c]
            cell = vtk_cells[offset + 1 + c : offset + 1 + c + numnodes]
            c += 1
            cell = (
                cell
                if cell_type not in pixel_voxel
                else cell[[0, 1, 3, 2]]  # type: ignore[call-overload]
                if cell_type == 8
                else cell[[0, 1, 3, 2, 4, 5, 7, 6]]  # type: ignore[call-overload]
            )
            cell_type = cell_type if cell_type not in pixel_voxel else cell_type + 1
            cell_type = vtk_to_meshio_type[cell_type]

        if len(cells) > 0 and cells[-1][0] == cell_type:
            cells[-1][1].append(cell)
        else:
            cells.append((cell_type, [cell]))
    return cells
        


def make_tubes(model:str):
    config = read_config(f"configfiles/{model}.yml")
    dt,T = config["dt"], config["T"]
    times = np.arange(0, T + dt, dt*config["output_frequency"])
    for netw in ["artery", "vein"]:
        tubes = get_tubes(get_result(model, "artery", times).ctp()).triangulate()
        #pv.save_meshio(f"{netw}_tubes.xdmf", tubes)
        with meshio.xdmf.TimeSeriesWriter(f"{netw}_tubes.xdmf") as writer:
            writer.write_points_cells(tubes.points, 
                                      get_cells(tubes))
            for t in times:
                writer.write_data(t, point_data={"c": tubes[f"c_{t}"]})


if __name__ == "__main__":
    typer.run(make_tubes)
