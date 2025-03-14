import dolfin
import pyvista
from typing import Tuple, Union, List
import numpy.typing as npt
import numpy as np
import pathlib
from typing import Optional

# taken from https://github.com/jorgensd/dolfin_pyvista_adapter

_vtk_perm = {
    "interval": {i: np.arange(i+1) for i in range(1, 10)},
    "vertex": {1: 1},
    "triangle": {
        1: [0, 1, 2],
        2: [0, 1, 2, 5, 3, 4],
        3: [0, 1, 2, 7, 8, 3, 4, 6, 5, 9],
        4: [0, 1, 2, 9, 10, 11, 3, 4, 5, 8, 7, 6, 12, 13, 14],
    },
    "tetrahedron": {
        1: [0, 1, 2, 3],
        2: [0, 1, 2, 3, 9, 6, 8, 7, 5, 4],
        3: [0, 1, 2, 3, 14, 15, 8, 9, 13, 12, 10, 11, 6, 7, 4, 5, 18, 16, 17, 19],
    },
}


def create_vtk_structures(
    Vh: dolfin.FunctionSpace,
) -> Tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64]]:
    """Given a (discontinuous) Lagrange space, create pyvista compatible structures based on the
    dof coordinates

    :param Vh: the function space
    :return: Mesh topology, cell types and mesh geometry arrays
    """
    family = Vh.ufl_element().family()
    mesh = Vh.mesh()
    if not (family in ["Discontinuous Lagrange", "Lagrange", "DQ", "Q", "DP", "P"]):
        raise RuntimeError(
            "Can only create meshes from continuous or discontinuous Lagrange spaces"
        )

    num_cells = mesh.num_cells()
    if num_cells == 0:
        return (
            np.zeros(0, dtype=np.int32),
            np.zeros(0, dtype=np.int32),
            np.zeros((0, 3), dtype=np.float64),
        )
    d_cell = mesh.ufl_cell().cellname()
    if d_cell == "vertex":
        cell_type = 1
    elif d_cell == "interval":
        cell_type = 68
    elif d_cell == "triangle":
        cell_type = 69
    elif d_cell == "tetrahedron":
        cell_type = 71
    else:
        raise RuntimeError(f"Unsupported {d_cell=}")
    cell_types = np.full(num_cells, cell_type, dtype=np.int32)
    bs = 1 if Vh.num_sub_spaces() == 0 else Vh.num_sub_spaces()
    Vh = Vh if bs == 1 else Vh.sub(0).collapse()
    num_dofs_per_cell = len(Vh.dofmap().cell_dofs(0))
    topology = np.zeros((num_cells, num_dofs_per_cell + 1), dtype=np.int32)
    topology[:, 0] = num_dofs_per_cell
    dofmap = Vh.dofmap()
    num_vertices_local = dofmap.index_map().size(dofmap.index_map().MapSize.ALL)

    x = np.zeros((num_vertices_local, 3), dtype=np.float64)
    el = Vh.element()
    for i in range(num_cells):
        topology[i, 1:] = dofmap.cell_dofs(i)
        cell = dolfin.Cell(mesh, i)
        cell_dof_coords = el.tabulate_dof_coordinates(cell)
        for j, d in enumerate(dofmap.cell_dofs(i)):
            x[d, : cell_dof_coords.shape[1]] = cell_dof_coords[j]

    degree = Vh.ufl_element().degree()
    try:
        perm = _vtk_perm[d_cell][degree]
    except KeyError:
        raise RuntimeError(
            f"Unsupported plotting of space {family} of {degree=} on {d_cell}"
        )
    topology[:, 1:] = topology[:, 1:][:, perm]
    return topology, cell_types, x


def plot(
    uh: dolfin.Function,
    filename: Optional[pathlib.Path] = None,
    show_edges: bool = True,
    glyph_factor: float = 1,
    off_screen: bool = True,
    cmap: Union[str, List[str]] = "viridis",
):
    """
    Plot a (discontinuous) Lagrange function with Pyvista

    :param uh: The function
    :param filename: If set, writes the plot to file instead of displaying it interactively
    :param show_edges: Show mesh edges if ``True``
    :param glyph_factor: Scaling of glyphs if input function is a function from a
      ``dolfin.VectorFunctionSpace``.
    :param off_screen: If ``True`` generate plots with virtual frame buffer using ``xvfb``.
    :param cmap: The colormap to use. If a string uses matplotlib colormaps,
      if a list uses strings as colors.
    """
    Vh = uh.function_space()

    dofmap = Vh.dofmap()
    num_vertices_local = dofmap.index_map().size(dofmap.index_map().MapSize.ALL)
    bs = 1 if Vh.num_sub_spaces() == 0 else Vh.num_sub_spaces()

    u_vec = uh.vector().get_local(np.arange(num_vertices_local))
    if bs > 1:
        u_out = np.zeros((len(u_vec) // bs, 3), dtype=np.float64)
        u_out[:, :bs] = u_vec.reshape(len(u_vec) // bs, bs)
    else:
        u_out = u_vec
    topology, cell_types, x = create_vtk_structures(Vh)

    grid = pyvista.UnstructuredGrid(topology, cell_types, x)
    grid[uh.name()] = u_out

    if bs > 1:
        grid.set_active_scalars(None)

    pyvista.OFF_SCREEN = off_screen
    pyvista.start_xvfb()

    plotter = pyvista.Plotter()

    # Workaround for pyvista issue for higherorder grids
    degree = Vh.ufl_element().degree()
    if degree > 1 and show_edges:
        first_order_el = Vh.ufl_element().reconstruct(degree=1)
        P1 = dolfin.FunctionSpace(Vh.mesh(), first_order_el)
        t1, c1, x1 = create_vtk_structures(P1)
        p1_grid = pyvista.UnstructuredGrid(t1, c1, x1)
        if bs > 1:
            plotter.add_mesh(p1_grid, show_edges=show_edges, cmap=cmap)
        else:
            plotter.add_mesh(grid, cmap=cmap)
            plotter.add_mesh(p1_grid, style="wireframe", cmap=cmap)
    else:
        plotter.add_mesh(grid, show_edges=show_edges, cmap=cmap)

    if bs > 1:
        glyphs = grid.glyph(orient=uh.name(), factor=glyph_factor)
        plotter.add_mesh(glyphs, cmap=cmap)
    print(filename)
    if filename is None:
        plotter.show()
    else:
        plotter.screenshot(filename)