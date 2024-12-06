from pubfig import (Units, FigureSpec, PanelsSpec, ElemSize, Panel, PanelFig, VectorImage, 
                    Location, RasterImage, compositor, AutoLabelOptions, Text, ImageType)
from pathlib import Path

units = Units.cm
fs = 4.9
fac = 0.9
xstart = 0.1
x2 = 10.4
x2r2 = 8.2
x2w = 3.4
textargs = {"x":1, "y":0.4, "size":8, "anchor":"left",
             "weight":"normal", "font":"DejaVu Sans"}
textargs_bar = dict(textargs, y=0)
offset = Location(-0.4,-0.3)
y0 = 0.4
y1 = 8
y2 = 11.4
y3 = 16.2
y4 = 21.2

label_bar_size=(4.3,4)
vel_hist_size=(4,3.6)

class Figure2_Panels(PanelsSpec):
    a: Panel = Panel(
            RasterImage("plots/meshplots/labeled_arteries.png", ElemSize(10.8, 7, units)),
            Location(xstart, y0),
            content_offset=Location(-0.4,-0.0),
            text=Text("Arterial network", **dict(textargs, x=4.5, y=0))
        )
    b: Panel = Panel(
            RasterImage("paper/figures/peristaltic_flow.png", ElemSize(7, 2.7, units)),
            Location(x2, y0),
            text=Text("Peristaltic flow", **dict(textargs, x=2, y=0)),
            content_offset=Location(0.5,-0.9),
        )
    c: Panel = Panel(
        RasterImage("plots/pvs_flow_prod/sas_flow-arteries/sas_flow-arteries_radius.png",
                     ElemSize(*label_bar_size, units)),
        Location(x2, 2.0),
        content_offset=Location(-0.5,-0.5),
        text=Text("radii of selected arteries", **dict(textargs, y=0.55, x=1))
    )
    c2: Panel = Panel(
        RasterImage("plots/pvs_flow_prod/sas_flow-arteries/sas_flow-arteries_velocity.png",
                     ElemSize(*label_bar_size, units)),
        Location(x2, 5),
        content_offset=Location(-0.5,-0.5),
        text=Text("pressure-driven PVS flow", **dict(textargs, y=0.55, x=1))
    )
    c3: Panel = Panel(
        RasterImage("plots/pvs_flow_peristaltic/cardiac_pvs_oscillation/cardiac_pvs_oscillation_velocity.png",
                     ElemSize(*label_bar_size, units)),
        Location(x2  + 4, 2.0),
        content_offset=Location(-0.5,-0.5),
        text=Text("cardiac-driven PVS flow", **dict(textargs, y=0.55, x=1))
    )
    c4: Panel = Panel(
        RasterImage("plots/pvs_flow_peristaltic/vasomotion-strong/vasomotion-strong_velocity.png",
                     ElemSize(*label_bar_size, units)),
        Location(x2 + 4, 5),
        content_offset=Location(-0.5,-0.5),
        text=Text("vasomotion-driven PVS flow", **dict(textargs, y=0.55, x=1))
    )
    m: Panel = Panel(
        RasterImage("plots/modelA/modelA_conc_at_label.png",
                     ElemSize(8, 3, units)),
        Location(xstart, y1),
        #content_offset=Location(-1,-1.3),
        text=Text("model A", **dict(textargs, y=0.5))
    )
    n: Panel = Panel(
        RasterImage("plots/modelA-strongVM/modelA-strongVM_conc_at_label.png",
                     ElemSize(8, 4, units)),
        Location(xstart, y1 + 2.7), auto_label=False,
        #content_offset=Location(-1,-1.3),
        text=Text("model A + VM", **dict(textargs, y=0.5))
    )
    d: Panel = Panel(
        RasterImage("plots/pvs_flow_prod/sas_flow-arteries/sas_flow-arteries_velocity_histo_cell.png",
                     ElemSize(*vel_hist_size, units)),
        Location(x2r2, y1),
        content_offset=Location(-0.3,0),
        text=Text("pressure-driven", **dict(textargs, y=0.1))
    )
    e: Panel = Panel(
        RasterImage("plots/pvs_flow_peristaltic/cardiac_pvs_oscillation/cardiac_pvs_oscillation_velocity_histo_cell.png",
                     ElemSize(*vel_hist_size, units)),
        Location(x2r2 + x2w, y1), auto_label=False,
        content_offset=Location(-0.3,0),
        text=Text("cardiac-driven", **dict(textargs, y=0.1))
    )
    f: Panel = Panel(
        RasterImage("plots/pvs_flow_peristaltic/vasomotion-strong/vasomotion-strong_velocity_histo_cell.png",
                     ElemSize(*vel_hist_size, units)),
        Location(x2r2 + 2*x2w, y1), auto_label=False,
        content_offset=Location(-0.3,0),
        text=Text("vasomotion-driven", **dict(textargs, y=0.1))
    )
    i: Panel = Panel(
        RasterImage("plots/modelA/modelA_ridgeline_total_smoothed.png",
                     ElemSize(4, 4.5, units)),
        Location(x2r2, y2),
        #content_offset=Location(-1,-1.3),
        #text=Text("model A total PVS tracer content", **dict(textargs, y=-0.3))
    )
    j: Panel = Panel(
        RasterImage("plots/modelA-PVS-disp/modelA-PVS-disp_ridgeline_total_smoothed.png",
                     ElemSize(4, 4.5, units)),
        Location(x2r2 + x2w, y2),
        #content_offset=Location(-1,-1.3),
        #text=Text("model A + PVS disp (10x) - total PVS tracer content", **dict(textargs, y=-0.3))
    )
    j2: Panel = Panel(
        RasterImage("plots/modelA-strongVM/modelA-strongVM_ridgeline_total_smoothed.png",
                     ElemSize(4, 4.5, units)),
        Location(x2r2 + 2*x2w, y2),
        #content_offset=Location(-1,-1.3),
        #text=Text("model A + VM total PVS tracer content", **dict(textargs, y=-0.3))
    )
    m1: Panel = Panel(
        RasterImage("plots/modelA/modelA_1-2-3-4_MCA-R-0.1-2.0-8.0-12.0_details.png",
                     ElemSize(8, 2.0, units),crop=(0,0,0,150)),
        Location(xstart, y1+6.0),
        content_offset=Location(0,-0.1),
        #text=Text("model A", **dict(textargs, y=0.5))
    )
    m2: Panel = Panel(
        RasterImage("plots/modelA-strongVM/modelA-strongVM_1-2-3-4_MCA-R-0.1-2.0-8.0-12.0_details.png",
                     ElemSize(8, 2.2, units), crop=(0,100,0,0)),
        Location(xstart, y1 +6.0 + 1.8), auto_label=False,
        content_offset=Location(0,-0.2),
        #text=Text("model A + VM", **dict(textargs, y=0.5))
    )
    m3: Panel = Panel(
        RasterImage("plots/modelA/modelA_2-4-6-8_ACA-A4-0.1-0.2-1.0-2.0_details.png",
                     ElemSize(8, 2.0, units),crop=(0,0,0,150)),
        Location(xstart, y1+9.6),
        content_offset=Location(0,-0.1),
        #text=Text("model A", **dict(textargs, y=0.5))
    )
    m4: Panel = Panel(
        RasterImage("plots/modelA-strongVM/modelA-strongVM_2-4-6-8_ACA-A4-0.1-0.2-1.0-2.0_details.png",
                     ElemSize(8, 2.2, units), crop=(0,100,0,0)),
        Location(xstart, y1 +9.6 + 1.8), auto_label=False,
        content_offset=Location(0,-0.2),
        #text=Text("model A + VM", **dict(textargs, y=0.5))
    )

    l: Panel = Panel(
        RasterImage("plots/comparisons/modelA_modelA-strongVM_modelA-PVS-disp/modelA_modelA-strongVM_modelA-PVS-disp_fta.png",
                     ElemSize(10, 5, units)),
        Location(x2r2, y3),
        #content_offset=Location(-1,-1.3),
        #text=Text("model A + VM total PVS tracer content", **dict(textargs, y=-0.3))
    )
    k: Panel = Panel(
        RasterImage("plots/modelA/modelA_overview_4-6.png",
                     ElemSize(6, 8.5, units), crop=(50, 0, 630, 0)),
        Location(xstart, y4),
        content_offset=Location(-0.3, 0),
        text=Text("model A", **dict(textargs, y=0.2, x= 3))
    )
    k1: Panel = Panel(
        RasterImage("plots/modelA-PVS-disp/modelA-PVS-disp_overview_4-6.png",
                     ElemSize(6, 8.5, units), crop=(50, 0, 630, 0)),
        Location(xstart + 4.9, y4),auto_label=False,
        #content_offset=Location(-1,-1.3),
        text=Text("model A + PVS disp", **dict(textargs, y=0.2, x= 3))
    )
    k2: Panel = Panel(
        RasterImage("plots/modelA-strongVM/modelA-strongVM_overview_4-6.png",
                     ElemSize(8, 8.5, units), crop=(50, 0, 20, 0)),
        Location(xstart + 10.3, y4), auto_label=False,
        #content_offset=Location(-1,-1.3),
        text=Text("model A + VM", **dict(textargs, y=0.2, x= 3))
    )

class Figure2(FigureSpec):
    figure_size = ElemSize(18.5, 29.8, units)
    output_file = Path.cwd() / "plots" / "figures" / "figure3.pdf"
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