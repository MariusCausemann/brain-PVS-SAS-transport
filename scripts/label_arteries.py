import pyvista as pv
import numpy as np
import seaborn as sns
from cmap import Colormap
import pandas as pd
from PIL import Image, ImageFont, ImageDraw  
from test_map_on_global_coords_shift import map_kdtree

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
   art = pv.read("mesh/networks/arteries_tubes.vtk").triangulate()

   art["radius"] *= 1e3
   col = "cyan"   
   cmap = Colormap("curl_pink_r")   
   arrow_length = 0.015
   window_size = (1200, 850)
   camera_direction = (-0.256, -0.958, -0.122)
   #pointlabels.sort(key=lambda pl: pl[1][2], reverse=True)
   #from IPython import embed; embed()
   idx = map_kdtree(np.array(art.points), 
                    np.array([p for l,p in pointlabels]).astype(np.float32),
                     distance_upper_bound=3e-3, k=100)

   ex = art.extract_points(np.unique(idx))
   #from IPython import embed; embed()
   bar_args=dict(title="radius (mm)", vertical=False,
                     height=0.07, width=0.6, position_x=0.2,
                     position_y=0.0, title_font_size=28,
                     label_font_size=28, fmt="%.1f")
   pl = pv.Plotter(off_screen=True, window_size=window_size)
   pl.add_mesh(art, cmap=cmap, scalar_bar_args=bar_args)
   pl.add_mesh(ex, show_scalar_bar=False, color=col)
   pl.camera_position = 'xz'
   #pl.camera.roll += 5
   pl.camera.azimuth += 165
   pl.camera.elevation += 10
   pl.camera.zoom(1.85)
   pl.camera.focal_point = np.array(pl.camera.focal_point) + np.array([0,0,-0.008])

   import vtk
   coordinate = vtk.vtkCoordinate()
   coordinate.SetCoordinateSystemToWorld()

   df = pd.DataFrame(pointlabels, columns=["label", "coord"])

   def get_disp_coord(p):
      coordinate.SetValue(*p)
      return coordinate.GetComputedViewportValue(pl.renderer)
   df["xy"] = df["coord"].apply(get_disp_coord)
   df["left"] = df["xy"].apply(lambda xy: xy[0] > 620)
   df["y"] = df["xy"].apply(lambda xy: xy[1])

   def add_labels(df, pos, offset=0): 
      y = 0
      for _, row in df.sort_values("y").iterrows():
         y = max([y + 50, row["y"]])
         act = pl.add_text(row["label"], position=[pos[0], y], font_size=16)
         c_world = act.GetActualPositionCoordinate().GetComputedWorldValue(pl.renderer)
         #l = pv.Line(c_world, row["coord"])
         vec = -np.array(c_world) + np.array(row["coord"])
         l = pv.Arrow(c_world, direction=vec,
                     scale=np.linalg.norm(vec)*0.98, shaft_radius=2e-3, tip_radius=0.8e-2,
                     tip_length=5e-2)
         pl.add_mesh(l, color="black")

   #for _, row in df.iterrows():
   #   pl.add_points(np.array(row["coord"]), render_points_as_spheres=True,
   #                  point_size=20, opacity=0.5, color="yellow")

   #add_labels(df[df["left"]], pos=np.array([1020, 200]))
   #add_labels(df[df["left"]==False], pos=np.array([20, 200]))

   def add_labels2d(ax, df, xpos=0, cs=None, offset=0):
      y = 0 if offset>0 else df["y"].max()
      op = max if offset>0 else min
      for _, row in df.sort_values("y", ascending=True if offset>0 else False).iterrows():
         y = op([y + offset, row["y"] + offset])
         ax.annotate(text=row["label"],xy=row["xy"], 
                     arrowprops=dict(facecolor='black',#width=1, shrink=0.02,headwidth=4,headlength=4,
                                     arrowstyle="->",linewidth=2,
                                       connectionstyle=cs),
                      xytext=(xpos, y), fontsize=16)
         

   img = pl.screenshot(filename, transparent_background=True)
   import matplotlib.pyplot as plt
   fig, ax = plt.subplots(figsize=(16,12))
   #fig.figimage(img, resize=True)
   ax.set_axis_off()
   ax.imshow(img, origin='upper', extent=[0, window_size[0], 0, window_size[1]])
   add_labels2d(ax, df[(df["left"]) & (df["y"] > 300)], xpos=1100, cs ="angle,angleA=0,angleB=-110,rad=50", offset=40)
   add_labels2d(ax, df[(df["left"]) & (df["y"] < 300)], xpos=1100, cs ="angle,angleA=0,angleB=110,rad=50",offset=-40)

   add_labels2d(ax, df[(df["left"]==False) & (df["y"] > 120)], xpos=20, cs ="angle,angleA=0,angleB=-70,rad=50", offset=40)
   add_labels2d(ax, df[(df["left"]==False) & (df["y"] < 120)], xpos=20, cs ="angle,angleA=0,angleB=70,rad=50", offset=-40)

   ax.set_xmargin(0);ax.set_ymargin(0)
   plt.savefig(filename, transparent=True, pad_inches=0.0,bbox_inches='tight', dpi=300)

   print(pl.camera.direction)
   pl.export_html(filename.replace("png", "html"))
