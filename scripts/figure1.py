from pubfig import (Units, FigureSpec, PanelsSpec, ElemSize, Panel, PanelFig, VectorImage, 
                    Location, RasterImage, compositor, AutoLabelOptions, Text, ImageType)
from pathlib import Path

units = Units.cm
xstart = 0.1
textargs = {"x":2, "y":0.4, "size":8, "anchor":"middle",
             "weight":"normal", "font":"DejaVu Sans"}
textargs_histo = dict(textargs, y=0)
y0 = 0.4
y1 = 7.4
w = 2.65


class Figure2_Panels(PanelsSpec):
    a: Panel = Panel(
            RasterImage("paper/figures/Brain-PVS-callouts.png",
                         ElemSize(9, 7, units),),
            Location(xstart, y0),
            content_offset=Location(0,-0.5),
        )
    b: Panel = Panel(
            RasterImage("plots/modelA/modelA_overview.png",
                         ElemSize(13, 8, units), crop=(0,100,0,0)),
            Location(xstart, y1),
            content_offset=Location(0,0),
            text=(Text(" 1:00 h", **dict(textargs, y=-0.1, x=w/2)),
                  Text(" 6:00 h", **dict(textargs, y=-0.1, x=3*w/2 )),
                  Text("12:00 h", **dict(textargs, y=-0.1, x=5*w/2)),
                  Text("24:00 h", **dict(textargs, y=-0.1, x=7*w/2)))
        )


class Figure2(FigureSpec):
    figure_size = ElemSize(18.2, 16, units)
    output_file = Path.cwd() / "plots" / "figures" / "figure1.png"
    auto_label_options = AutoLabelOptions(
        first_char:=Text("A", 0, 0, size=12, weight="bold"))
    panels = Figure2_Panels()
    generate_image = ImageType.png  # Could also be `png` or `none`
    image_dpi = 300
    #plot_grid_every = 1

@compositor(Figure2, memoize_panels=True, recompute_panels=False)
def create_fig1(figure: Figure2):
    pass

create_fig1(recompute_panels=False)