import pyvista as pv
import numpy as np
import seaborn as sns
pointlabels = [("BA",(0.084971, 0.104599, 0.0424848)),
               ("ICA-R",(0.0965484, 0.112456,0.0413785)),
               ("ICA-L",(0.0695824, 0.111343, 0.0415065)),
               ("MCA-R",(0.0975131, 0.110452, 0.0626856)),
               ("MCA-L",(0.063242, 0.115508, 0.0613588)),
               ("MCA2-R",(0.127254, 0.118198, 0.0756301)),
               ("MCA2-L",(0.0360319, 0.114267, 0.0717319)),
               ("PCA-R",(0.100641, 0.0843866, 0.0609018)),
               ("PCA-L",(0.0686598, 0.0831056, 0.0568451)),
               ("ACA-A1-R", (0.0864598, 0.112447, 0.0623659)),
               ("ACA-A1-L", (0.0746942, 0.109032, 0.0643898)),
               ("ACA-A2", (0.0805165,0.115209,0.0781812)),
               ("ACA-A3", (0.0795565, 0.133436, 0.0905888)),
               ("ACA-A4", (0.0815724, 0.111484, 0.107046)),
               #("ACA-A5", (0.0815178, 0.0984657, 0.167766)),
               ("PER-R", (0.123073, 0.101526, 0.135256)),
               ("PER-L", (0.0401122, 0.0939796, 0.13324))
            ]

def get_tangent(v, n):
    return v - np.dot(np.dot(v,n), n)

if __name__=="__main__":
   filename = "plots/meshplots/labeled_arteries.png"
   art = pv.read("mesh/networks/arteries_tubes.vtk")
   labelcoords = lambda label: art.bounds[::2]

   art["radius"] *= 1e3
   col = "fuchsia"      
   arrow_length = 0.015
   camera_direction = (-0.256, -0.958, -0.122)
   points = [p for l, p in pointlabels]
   labels = [l for l, p in pointlabels]
   bar_args=dict(title="radius (mm)", vertical=False,
                     height=0.07, width=0.6, position_x=0.2,
                     position_y=0.0, title_font_size=24,
                     label_font_size=24, fmt="%.1f")
   pl = pv.Plotter(off_screen=True)
   pl.add_mesh(art, cmap="coolwarm", scalar_bar_args=bar_args, show_scalar_bar=True)
   pl.add_points(np.array(points), render_points_as_spheres=False, 
                 point_size=20, color=col)
   center = art.center_of_mass() - np.array([0,0, 0.047])
   label_coords = []
   for i, (l, p) in enumerate(pointlabels):
      dir = np.array(p) - center
      dir = get_tangent(dir, camera_direction)
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
   pl.screenshot(filename, transparent_background=True)
   print(pl.camera.direction)
   pl.export_html(filename.replace("png", "html"))
