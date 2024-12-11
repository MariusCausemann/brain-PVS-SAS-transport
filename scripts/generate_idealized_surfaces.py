import os 
import pyvista as pv
import numpy as np

os.makedirs("mesh/idealized/surfaces/", exist_ok=True)
r = 0.1
z = np.array([0,0,1])
skull = pv.Box((0.0103891, 0.155952,
                 0.0114499, 0.173221, 
                 0.001, 0.154949))
c = skull.center_of_mass()
par = pv.Sphere(r*0.6, center=c)
LV = pv.Sphere(r*0.2, center=c)
V34 = pv.Cylinder(c -0.4*r*z, direction=z, height=r*0.42,
                   radius=0.06*r)
ventricles = pv.merge([LV, V34])

for s,n in zip([V34, LV, par, skull, ventricles],
                ["V34", "LV", "parenchyma_incl_ventr",
                  "skull", "ventricles"]):
    pv.save_meshio(f"mesh/idealized/surfaces/{n}.ply",s)

# unused, but maybe useful one day
l = pv.lines_from_points([c -r*z, c -0.5*r*z, c])
l["radius"] = 0.01*np.ones(l.n_points)
l["root"] = [2, 0, 1]
l.save("mesh/networks/simple_line.vtk")