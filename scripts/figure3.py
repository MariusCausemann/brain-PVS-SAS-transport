from pubfig import (Units, FigureSpec, PanelsSpec, ElemSize, Panel, PanelFig, VectorImage, 
                    Location, RasterImage, compositor, AutoLabelOptions, Text, ImageType)
from pathlib import Path

units = Units.cm
fs = 4.9
fac = 0.9
xstart = 0.1
textargs = {"x":fs/2, "y":0.4, "size":8, "anchor":"middle",
             "weight":"normal", "font":"DejaVu Sans"}
textargs_histo = dict(textargs, y=0)
offset = Location(-0.4,-0.3)
y0 = 0.4
y1 = 7
y2 = 12
y3 = 18
y4 = 24
y5 = 33
y6 = 40
y7 = 44

class Figure2_Panels(PanelsSpec):
    a: Panel = Panel(
            RasterImage("plots/meshplots/labeled_arteries.png", ElemSize(10.8, 7, units)),
            Location(xstart, y0),
            content_offset=Location(0,-0.5),
            text=Text("Arterial network", **dict(textargs, x=4.5, y=0))
        )
    a2: Panel = Panel(
            RasterImage("paper/figures/peristaltic_flow.png", ElemSize(5, 2, units)),
            Location(xstart + 11, y0),
            text=Text("Peristaltic flow", **dict(textargs, x=2, y=0))
        )
    b: Panel = Panel(
        RasterImage("plots/pvs_flow_prod/sas_flow-arteries/sas_flow-arteries_radius.png",
                     ElemSize(6, 6, units)),
        Location(xstart + 11, 2.5),
        content_offset=Location(-0.5,-0.8),
        text=Text("radii of selected arteries", **dict(textargs, y=1.2, x=3))
    )
    c: Panel = Panel(
        RasterImage("plots/pvs_flow_prod/sas_flow-arteries/sas_flow-arteries_velocity_histo_cell.png",
                     ElemSize(6, 6, units)),
        Location(xstart, y1),
        #content_offset=Location(-1,-1.3),
        text=Text("pressure-driven PVS flow velocity", **dict(textargs, y=0.3))
    )
    d: Panel = Panel(
        RasterImage("plots/pvs_flow_peristaltic/cardiac_pvs_oscillation/cardiac_pvs_oscillation_velocity_histo_cell.png",
                     ElemSize(6, 6, units)),
        Location(xstart + 6, y1),
        #content_offset=Location(-1,-1.3),
        text=Text("cardiac peristaltic PVS flow velocity", **dict(textargs, y=0.3))
    )
    e: Panel = Panel(
        RasterImage("plots/pvs_flow_peristaltic/vasomotion-strong/vasomotion-strong_velocity_histo_cell.png",
                     ElemSize(6, 6, units)),
        Location(xstart + 12, y1),
        #content_offset=Location(-1,-1.3),
        text=Text("vasomotion peristaltic PVS flow velocity", **dict(textargs, y=0.3))
    )
    f: Panel = Panel(
        RasterImage("plots/pvs_flow_prod/sas_flow-arteries/sas_flow-arteries_velocity.png",
                     ElemSize(6, 6, units)),
        Location(xstart, y2),
        #content_offset=Location(-1,-1.3),
        text=Text("pressure-driven PVS flow velocity", **dict(textargs, y=0.3))
    )
    g: Panel = Panel(
        RasterImage("plots/pvs_flow_peristaltic/cardiac_pvs_oscillation/cardiac_pvs_oscillation_velocity.png",
                     ElemSize(6, 6, units)),
        Location(xstart + 6, y2),
        #content_offset=Location(-1,-1.3),
        text=Text("cardiac peristaltic PVS flow velocity", **dict(textargs, y=0.3))
    )
    h: Panel = Panel(
        RasterImage("plots/pvs_flow_peristaltic/vasomotion-strong/vasomotion-strong_velocity.png",
                     ElemSize(6, 6, units)),
        Location(xstart + 12, y2),
        #content_offset=Location(-1,-1.3),
        text=Text("vasomotion peristaltic PVS flow velocity", **dict(textargs, y=0.3))
    )
    i: Panel = Panel(
        RasterImage("plots/modelA/modelA_ridgeline_total_smoothed.png",
                     ElemSize(6, 6, units)),
        Location(xstart, y3),
        #content_offset=Location(-1,-1.3),
        text=Text("model A total PVS tracer content", **dict(textargs, y=-0.3))
    )
    j: Panel = Panel(
        RasterImage("plots/modelA-PVS-disp/modelA-PVS-disp_ridgeline_total_smoothed.png",
                     ElemSize(6, 6, units)),
        Location(xstart + 6, y3),
        #content_offset=Location(-1,-1.3),
        text=Text("model A + PVS disp (10x) - total PVS tracer content", **dict(textargs, y=-0.3))
    )
    j2: Panel = Panel(
        RasterImage("plots/modelA-strongVM/modelA-strongVM_ridgeline_total_smoothed.png",
                     ElemSize(6, 6, units)),
        Location(xstart + 12, y3),
        #content_offset=Location(-1,-1.3),
        text=Text("model A + VM total PVS tracer content", **dict(textargs, y=-0.3))
    )
    k: Panel = Panel(
        RasterImage("plots/modelA-strongVM/modelA-strongVM_overview.png",
                     ElemSize(14, 9, units)),
        Location(xstart, y4),
        #content_offset=Location(-1,-1.3),
        text=Text("model A + VM", **dict(textargs, y=0.2, x= 3))
    )
    l: Panel = Panel(
        RasterImage("plots/comparisons/modelA_modelA-strongVM_modelA-PVS-disp/modelA_modelA-strongVM_modelA-PVS-disp_fta.png",
                     ElemSize(14, 7, units)),
        Location(xstart, y5),
        #content_offset=Location(-1,-1.3),
        #text=Text("model A + VM total PVS tracer content", **dict(textargs, y=-0.3))
    )
    m: Panel = Panel(
        RasterImage("plots/modelA/modelA_conc_at_label.png",
                     ElemSize(14, 5, units)),
        Location(xstart, y6),
        #content_offset=Location(-1,-1.3),
        text=Text("model A", **dict(textargs, y=0.5))
    )
    nm: Panel = Panel(
        RasterImage("plots/modelA-strongVM/modelA-strongVM_conc_at_label.png",
                     ElemSize(14, 5, units)),
        Location(xstart, y7), auto_label=False,
        #content_offset=Location(-1,-1.3),
        text=Text("model A + VM", **dict(textargs, y=0.5))
    )




class Figure2(FigureSpec):
    figure_size = ElemSize(18.2, 50, units)
    output_file = Path.cwd() / "plots" / "figures" / "figure3.png"
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