import pyvista as pv

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
   art["radius"] *= 1e3

   points = [p for l, p in pointlabels]
   labels = [l for l, p in pointlabels]
   bar_args=dict(title="radius (mm)", vertical=False,
                     height=0.07, width=0.6, position_x=0.2,
                     position_y=0.0, title_font_size=24,
                     label_font_size=24, fmt="%.1f")
   pl = pv.Plotter(off_screen=True)
   pl.add_mesh(art, cmap="coolwarm", scalar_bar_args=bar_args)
   pl.add_point_labels(points, labels, font_size=24,
   margin=4, shadow=True, shape="rounded_rect",shape_opacity=0.7,
   render_points_as_spheres=True,
   shape_color="white",
   #justification_horizontal="center", 
   always_visible=True,
   point_color="red", point_size= 18, show_points=True)
   pl.camera_position = 'zx'
   pl.camera.roll += 90
   pl.camera.azimuth -= 25
   pl.camera.zoom(1.6)
   pl.screenshot(filename, transparent_background=False)
