#import "utils.typ": setup-figure, panel-label, panel-title, inline-header, cropped-image, y-label
#show: setup-figure


// ==========================================
// MAIN FIGURE LAYOUT
// ==========================================
#figure(
  grid(
    columns: 1,
    row-gutter: 10pt, 

    // ------------------------------------------
    // TOP ROW: Panels A, B, and (C & D)
    // ------------------------------------------
    grid(
      columns: (1fr, 1fr, 2.5fr), // A is normal, B is narrow, C/D are wide
      column-gutter: 5pt,
      
      // Column 1: Panel A
      [
        #panel-label("A")
        #v(5pt)
        #panel-title("PVS dilation")
        #image("../paper/figures/dilated_pvs.png", width: 110%)
      ],
      
      // Column 2: Panel B (Stacked Bar Charts)
      [
        #panel-label("B")
        #v(-20pt)
        #grid(
          columns: 1, row-gutter: 3pt,
          [
            #panel-title("max PVS velocity")
            #image("../plots/comparisons/pvs_flow_peristaltic+cardiac_pvs_oscillation_pvs_flow_peristaltic+cardiac_pvs_oscillation_enlarged_pvs_flow_peristaltic+vasomotion-strong_pvs_flow_peristaltic+vasomotion-strong-enlarged/pvsumax_comp.png", width: 100%)
          ],
          [
            #panel-title("mean PVS velocity")
            #image("../plots/comparisons/pvs_flow_peristaltic+cardiac_pvs_oscillation_pvs_flow_peristaltic+cardiac_pvs_oscillation_enlarged_pvs_flow_peristaltic+vasomotion-strong_pvs_flow_peristaltic+vasomotion-strong-enlarged/pvsumean_comp.png", width: 100%)
          ]
        )
      ],
      
      // Column 3: Panels C and D (Side-by-side Brains)
      grid(
        columns: (1fr, 1fr),
        column-gutter: 0pt,
        [
          #panel-label("C")
          #v(-5pt)
          #panel-title("vasomotion induced PVS velocity (control)")
          #image("../plots/pvs_flow_peristaltic/vasomotion-strong/vel3D.png", width: 130%)
        ],
        [
          #panel-label("D")
          #v(-5pt)
          #panel-title("vasomotion induced PVS velocity (dilated)")
          #image("../plots/pvs_flow_peristaltic/vasomotion-strong-enlarged/vel3D.png", width: 130%)
        ]
      ),
    ),
    v(-30pt),
    // ------------------------------------------
    // MIDDLE ROW: Panels E and F
    // ------------------------------------------
    grid(
      columns: (1fr, 1.1fr), // Right side slightly wider for the line plot
      column-gutter: 0pt,
      
      // Column 1: Panel E (Stacked Details)
      [
        #panel-label("E")
        #v(-40pt),
        #grid(
          columns: (auto, 1fr), gutter: 5pt, align: horizon,
          
          y-label("high PVS flow", nudge: 0.49cm),
          cropped-image("../plots/modelA-strongVM/modelA-strongVM_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png", 2.3cm, 100%, 0cm, 0.55cm),
          y-label("high PVS flow -\ndilated", nudge: -0.22cm),
          cropped-image("../plots/LargePVS/LargePVS_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png", 2.3cm, 100%, 0cm, -0.22cm)
        )
      ],
      
      // Column 2: Panel F (Line Chart)
      [
        #panel-label("F")
        #v(-20pt)
        #panel-title("first-time arrival")
        #v(-10pt)
        #image("../plots/comparisons/modelA_LargePVSA_modelA-strongVM_LargePVS/modelA_LargePVSA_modelA-strongVM_LargePVS_fta.png", width: 100%)
      ]
    ),

    // ------------------------------------------
    // BOTTOM ROW: Panels G and H
    // ------------------------------------------
    v(-20pt),
    grid(
      columns: (1fr, 1.1fr), 
      column-gutter: 0pt,
      align: bottom, // Aligns the plots cleanly to the bottom baseline
      
      // Column 1: Panel G (Bar Chart)
      [
        #panel-label("G")
        #v(-35pt)
        #cropped-image("../plots/comparisons/modelA_LargePVSA_modelA-strongVM_LargePVS/modelA_LargePVSA_modelA-strongVM_LargePVS_tot_bar.png", 5.1cm, 100%, 0cm, -0.1cm, width: 100%)
      ],
      
      // Column 2: Panel H (Line Charts Grid)
      [
        #panel-label("H")
        #v(-15pt)
        #panel-title("mean tracer concentration")
        #v(-8pt)
        #cropped-image("../plots/comparisons/modelA_LargePVSA_modelA-strongVM_LargePVS/modelA_LargePVSA_modelA-strongVM_LargePVS_conc.png", 5.1cm, 100%, 0cm, 0.1cm, width: 110%)
      ]
    )
  )
)