
import numpy as np
from scipy.spatial import cKDTree
import pyvista as pv
from test_map_on_global_coords_shift import map_kdtree
from extract_vessels import pvnetwork_to_polydata, get_tubes
np.random.seed(0)
from tqdm import tqdm

def space_colonization(starting_points, attraction_points, influence_radius, 
                       growth_length,kill_dist):
    # Initialize the graph with the starting points
    points = starting_points.copy()
    lines = []

    #while len(attraction_points) > 0:
    for i in range(200):
        print("============")
        print(f" lines generated: {len(lines)}")
        print(f" remaining attraction points: {len(attraction_points)}")
        # Build a KD-tree for efficient nearest neighbor search
        tree = cKDTree(points)

        # Find the nearest tree point
        dist, idx = tree.query(attraction_points, distance_upper_bound=influence_radius)

        attr_point_sets = {tp:np.where(idx==tp) for tp in np.unique(idx)}
        if len(points) in attr_point_sets.keys():
            attr_point_sets.pop(len(points))

        for tp, aps in attr_point_sets.items():
            vec = np.mean(attraction_points[aps], axis=0) - points[tp]
            new_point = points[tp] + growth_length * vec / np.linalg.norm(vec)
            points.append(new_point)
            lines.append((tp, len(points) - 1))

        # Remove the connected attraction points
        tree = cKDTree(points)
        dist, idx = tree.query(attraction_points, distance_upper_bound=kill_dist, k=100)
        killed = []
        for i, indices in enumerate(idx):
            if not np.allclose(indices, len(points)):
                killed.append(i)
        attraction_points = np.delete(attraction_points, killed, axis=0)
        if len(attraction_points)==0: break

    return np.array(points), np.array(lines)


# Define the number of attraction points and the influence radius
num_attraction_points = 100000
influence_radius = 0.005
growth_length = 0.0002
kill_dist = 5*growth_length
murray_exponent = 3

cerebrum = pv.read("mesh/standard/surfaces/cerebrum.ply")
arteries = pv.read("mesh/networks/arteries_smooth.vtk")

attraction_points = np.random.uniform(low=cerebrum.bounds[::2], 
                                      high=cerebrum.bounds[1::2], 
                                      size=(num_attraction_points, 3)) 
attraction_points = pv.PolyData(attraction_points).compute_implicit_distance(cerebrum, inplace=True)
attraction_points = attraction_points.extract_points(attraction_points["implicit_distance"] < 0)

start_points = list(np.array(arteries.points[arteries["root"] == 1, :]))

# Call the space colonization function
points, lines = space_colonization(start_points, 
                                   attraction_points.points,
                                   influence_radius,
                                   kill_dist=kill_dist, growth_length=growth_length)

tree = pv.PolyData(points,lines=np.column_stack([2*np.ones(lines.shape[0], dtype=int), lines]))

orig_idx = map_kdtree(arteries.points, np.array(start_points).astype(np.float32))
ext_idx = map_kdtree(tree.points, np.array(start_points).astype(np.float32))

assert len(start_points) == len(orig_idx)
tree.point_data["radius"] = np.zeros(tree.n_points, dtype=np.float32)
tree["radius"][ext_idx] = arteries["radius"][orig_idx]

def extent_by_murray(netw, startidx, murray_exponent):
    for i in startidx:
        assert netw["radius"][i] > 0
        ns = [ni for ni in netw.point_neighbors(i) if netw["radius"][ni]==0]
        while len(ns) ==1:
            idx = ns.pop()
            netw["radius"][idx] = netw["radius"][i] 
            ns += [ni for ni in netw.point_neighbors(idx) if netw["radius"][ni]==0]

        netw["radius"][ns] = (netw["radius"][i]**murray_exponent / len(ns))**(1/murray_exponent)
        extent_by_murray(netw, ns, murray_exponent)

assert tree["radius"].sum() > 0 

branches = tree.split_bodies()

for i,br in tqdm(enumerate(branches)):
    extent_by_murray(br, np.where(br["radius"] > 0)[0], murray_exponent=murray_exponent)
    br.point_data["branch_id"] = i + 1 
tree = pv.merge(branches, merge_points=True)
assert (tree["radius"] > 0).all()
arteries.point_data["branch_id"] = 0
tree.save("tree.vtk")
merged = pv.merge([tree, arteries], merge_points=True)
tubes = get_tubes(merged)
tubes.save("tubes.vtk")
merged.save("merged.vtk")
attraction_points.save("points.vtk")