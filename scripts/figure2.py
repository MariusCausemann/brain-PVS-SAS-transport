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
y1 = 8.3
y2 = 13.4
y3 = 16
w = 2.65
hn = 8.0
compwidth = 2.5
hpx = 800
abfac = 0.9

class Figure2_Panels(PanelsSpec):
    a: Panel = Panel(
            RasterImage("plots/meshplots/subdomains.png", ElemSize(fs*abfac, fs*abfac, units)),
            Location(xstart, y0),
            content_offset=Location(0.2,-0.3),
            text=Text("3D model of SAS and ventricles", **textargs)
        )
    b: Panel = Panel(
            RasterImage("plots/meshplots/boundaries.png", ElemSize(fs*abfac, fs*abfac, units)),
            Location(xstart, y0 + fs*0.8),
            content_offset=Location(0.2,-0.3),
            text=Text("boundaries of the CSF space",  **textargs)
        )
    c: Panel = Panel(
            RasterImage("paper/figures/flow_and_dispersion.png", 
                        ElemSize(fs*1.65, fs*1.65, units)),
            Location(xstart + fs*0.95, y0),
            content_offset=Location(0.2,-0.3),
            text=Text("model illustration", **dict(textargs,x=fs/2*1.6))
        )
    d: Panel = Panel(
            RasterImage("results/csf_flow/sas_flow/csf_v.png", 
                        ElemSize(fs, fs, units)),
            Location(xstart + 3*fs*fac, y0),
            content_offset=offset,
            text=Text("production-driven CSF flow", **textargs)
    )
    d2: Panel = Panel(
            RasterImage("results/csf_flow/sas_flow/prod_velocity_histo.png",
                         ElemSize(fs, fs*0.6, units)),
            Location(xstart + 3*fs*fac, y0 + fs),
            content_offset=Location(-0.3, 0),
            text=Text("production-driven CSF velocities", **dict(textargs, y=0.1))
        )
    e: Panel = Panel(
            RasterImage("results/csf_flow/cardiac_sas_flow/csf_v.png", 
                        ElemSize(fs, fs, units)),
            Location(xstart, y1),
            content_offset=offset,
            text=Text("cardiac-driven CSF flow", **textargs)
        )
    g: Panel = Panel(
            RasterImage("results/csf_flow/cardiac_sas_flow/R.png",
                         ElemSize(fs, fs, units)),
            Location(xstart + fs*fac, y1),
            content_offset=offset,
            text=Text("cardiac dispersion enhancement", **textargs)
        )
    
    f: Panel = Panel(
            RasterImage("results/csf_flow/respiratory_sas_flow/csf_v.png", 
            ElemSize(fs, fs, units)),
            Location(xstart +fac*2*fs, y1),
            content_offset=offset,
            text=Text("respiratory-driven CSF flow", **textargs)
        )
    h: Panel = Panel(
            RasterImage("results/csf_flow/respiratory_sas_flow/R.png",
                         ElemSize(fs, fs, units)),
            Location(xstart +3*fs*fac, y1),
            content_offset=offset,
            text=Text("respiratory dispersion enhancement", **dict(textargs))
        )
    i: Panel = Panel(
            RasterImage("results/csf_flow/cardiac_sas_flow/cardiac_velocity_histo.png",
                         ElemSize(fs*fac, fs/2, units)),
            Location(xstart + 0*fs, y2),
            text=Text("cardiac flow velocities", **textargs_histo)
        )
    j: Panel = Panel(
            RasterImage("results/csf_flow/cardiac_sas_flow/R_histo.png",
                         ElemSize(fs*fac, fs/2, units)),
            Location(xstart + 1*fs*fac, y2),
            text=Text("cardiac dispersion enhancement", **textargs_histo)
        )
    k: Panel = Panel(
            RasterImage("results/csf_flow/respiratory_sas_flow/resp_velocity_histo.png",
                         ElemSize(fs*fac, fs/2, units)),
            Location(0.0 + 2*fs*fac, y2),
            text=Text("respiratory flow velocities", **textargs_histo)
        )
    l: Panel = Panel(
            RasterImage("results/csf_flow/respiratory_sas_flow/R_histo.png",
                         ElemSize(fs*fac, fs/2, units)),
            Location(0.0 +3*fs*fac, y2),
            text=Text("respiratory dispersion enhancement", **textargs_histo)
        )
    
    m: Panel = Panel(
            RasterImage("plots/modelA/modelA_overview.png",
                         ElemSize(12.9, hn, units), crop=(0,100,0,0)),
            Location(0.0, y3),
            content_offset=Location(0,0),
            text=(Text(" 1:00 h", **dict(textargs, y=-0.1, x=w/2)),
                  Text(" 6:00 h", **dict(textargs, y=-0.1, x=3*w/2 )),
                  Text("12:00 h", **dict(textargs, y=-0.1, x=5*w/2)),
                  Text("24:00 h", **dict(textargs, y=-0.1, x=7*w/2)))
        )
    n: Panel = Panel(
            RasterImage("plots/modelA-LowD/modelA-LowD_overview.png",
                         ElemSize(3, hn, units), crop=(1400, 100, 1300, 0)),
            Location(12.8, y3), auto_label=True,
            content_offset=Location(0,0),
            text=Text("low dispersion (12 h)", **dict(textargs, y=-0.1, x=1.5))
        )
    n2: Panel = Panel(
            RasterImage("plots/modelA-OnlyDispersion/modelA-OnlyDispersion_overview.png",
                         ElemSize(3, hn, units), crop=(1400, 100, 1300, 0)),
            Location(12.8 + compwidth, y3), auto_label=False,
            content_offset=Location(0,0),
            text=Text("no advection (12 h)", **dict(textargs, y=-0.1, x=1.5))
        )


class Figure2(FigureSpec):
    figure_size = ElemSize(18.2, 24, units)
    output_file = Path.cwd() / "plots" / "figures" / "figure2.png"
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