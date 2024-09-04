import pyvista as pv
filename = "plots/labeled_arteries.png"
art = pv.read("mesh/networks/arteries_tubes.vtk")
art["radius"] *= 1e3

pointlabels = [("ACA", (0.0836656391620636, 0.12674154341220856, 0.10104840993881226)),
               ("ICA",(0.09893926978111267, 0.11347689479589462, 0.03919529169797897)),
            ("ICA",(0.0687378594571459, 0.1128146693515864, 0.03797117973496493)),
           ("MCA",(0.11279711872339249, 0.12369046360254288, 0.0656440332531929)),
           ("MCA",(0.05993403121829033, 0.1203613206744194, 0.06074688211083412)),
            ("BA",(0.08326485485531022, 0.10635409315410013, 0.038923607855650535)),
            ("PCA",(0.06962447613477707, 0.057159677147865295, 0.07165350019931793)),
            ("PCA",(0.09846239537000656, 0.07383400201797485, 0.0722760483622551)),
            ]
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
