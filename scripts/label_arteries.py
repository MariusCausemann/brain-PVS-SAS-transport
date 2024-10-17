import pyvista as pv
import numpy as np
import seaborn as sns
pointlabels = [("ACA", (0.0797341, 0.132915, 0.0937325)),
               ("ICA1",(0.0942801, 0.11788, 0.0505869)),
               ("ICA2",(0.0679531, 0.118951, 0.0485859)),
               ("MCA1",(0.10513, 0.114464, 0.0672847)),
               ("MCA2",(0.0640877, 0.114662, 0.0613344)),
               ("BA",(0.0845697, 0.105035, 0.0450061)),
               ("PCA1",(0.0719072, 0.0640646, 0.063061)),
               ("PCA2",(0.0969292, 0.0740488, 0.0658269)),
            ]

if __name__=="__main__":
   filename = "plots/labeled_arteries.png"
   art = pv.read("mesh/networks/arteries_tubes.vtk")
   labelcoords = lambda label: art.bounds[::2]

   art["radius"] *= 1e3
   col = "fuchsia"      
   arrow_length = 0.03

   points = [p for l, p in pointlabels]
   labels = [l for l, p in pointlabels]
   bar_args=dict(title="radius (mm)", vertical=False,
                     height=0.07, width=0.6, position_x=0.2,
                     position_y=0.0, title_font_size=24,
                     label_font_size=24, fmt="%.1f")
   pl = pv.Plotter(off_screen=True)
   pl.add_mesh(art, cmap="coolwarm", scalar_bar_args=bar_args)
   pl.add_points(np.array(points), render_points_as_spheres=False, 
                 point_size=20, color=col)
   center = art.center_of_mass() - np.array([0,0, 0.047])
   label_coords = []
   for i, (l, p) in enumerate(pointlabels):
      dir = np.array(p) - center
      dir /= np.linalg.norm(dir)
      arr = pv.Arrow(start=np.array(p), direction=dir, scale=arrow_length,
                     tip_length=1e-3, tip_radius=5e-3, shaft_radius=1e-2)
      pl.add_mesh(arr, color=col)
      label_coords.append(p + arrow_length*dir)
      
   pl.add_point_labels(label_coords, labels, show_points=False, 
      always_visible=True, shape_color=col)
   pl.camera_position = 'xz'
   #pl.camera.roll += 10
   pl.camera.azimuth += 160
   pl.camera.elevation += 7
   pl.camera.zoom(1.6)
   pl.screenshot(filename, transparent_background=False)
