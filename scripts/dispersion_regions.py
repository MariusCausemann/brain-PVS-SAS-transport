import pyvista as pv
import numpy as np

grid = pv.read("mesh/T1/volmesh/mesh.xdmf")
grid = grid.extract_cells(np.isin(grid["label"], [1,3,4]))
n_regions = 100
idx = np.random.randint(0, grid.n_points, n_regions)
points = grid.points[idx]

#points = pv.Box(grid.bounds, level=4).points

def dist(grid, point):
    coords = grid.cell_centers().points
    d = grid.compute_implicit_distance(pv.Sphere(radius=1e-6,center=point)).ptc()
    return d["implicit_distance"]
    #return np.linalg.norm(coords - point, axis=1)

distances = np.array([dist(grid, p) for p in points])
grid["regions"] = np.argmin(distances, axis=0)
grid.pack_labels(scalars="regions", inplace=True)

grid.save("partioned.vtk")

p_tilde = 45.3 # Pa/m
x_tilde = 0.1
rho = 993 # kg/m^3
nu = 7e-7 # m^2/s
omega = 2*np.pi
h = 3e-3 / 2 
P = p_tilde / x_tilde /(rho*omega*nu/h)


D = 5.26e-10
Sc = nu / D
alpha = np.sqrt(h**2 * omega / nu)
beta2 = alpha**2 * Sc

R = P**2 / alpha**3