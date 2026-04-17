#import "utils.typ": setup-figure, panel-label, panel-title, inline-header, cropped-image, y-label, crop-image, x-label
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
      columns: (1fr, 1.5fr), // A is normal, B is narrow, C/D are wide
      column-gutter: 0pt,
      
      // Column 1: Panel A
      [
        #panel-label("A")
        #v(-20pt)
        #panel-title("PVS dilation")  #v(-6pt)
        #image("../paper/figures/dilated_pvs.png", width: 65%)
        #panel-label("B", nudge-dy: -10pt)
        #v(-20pt)
        #grid(
          columns: 2, row-gutter: 5pt,
          [
            #v(-2pt)
            #panel-title("Max PVS velocity")
            #v(-2pt)
            #image("../plots/comparisons/pvs_flow_peristaltic+cardiac_pvs_oscillation_pvs_flow_peristaltic+cardiac_pvs_oscillation_enlarged_pvs_flow_peristaltic+vasomotion-strong_pvs_flow_peristaltic+vasomotion-strong-enlarged/pvsumax_comp.png", width: 105%)
          ],
          [
            #v(-2pt)
            #panel-title("Mean PVS velocity")
            #v(-2pt)
            #image("../plots/comparisons/pvs_flow_peristaltic+cardiac_pvs_oscillation_pvs_flow_peristaltic+cardiac_pvs_oscillation_enlarged_pvs_flow_peristaltic+vasomotion-strong_pvs_flow_peristaltic+vasomotion-strong-enlarged/pvsumean_comp.png", width: 105%)
          ]
        )
      ],
      

      
      // Column 3: Panels C and D (Side-by-side Brains)
      grid(
        columns: (1fr, 1fr),
        column-gutter: 0pt,
        [
          #panel-label("C", nudge-dx: -5pt)
          #v(-16pt)
          #panel-title("Vasomotion induced PVS velocity (control)")
          #image("../plots/pvs_flow_peristaltic/vasomotion-strong/vel3D_black.png", width: 135%)
        ],
        [
          #panel-label("D", nudge-dx: -5pt)
          #v(-16pt)
          #panel-title("Vasomotion induced PVS velocity (dilated)")
          #image("../plots/pvs_flow_peristaltic/vasomotion-strong-enlarged/vel3D_black.png", width: 135%)
        ]
      ),
    ),
    v(-15pt),
    // ------------------------------------------
    // MIDDLE ROW: Panels E and F
    // ------------------------------------------
    grid(
      columns: (1fr, 1.1fr), // Right side slightly wider for the line plot
      column-gutter: 0pt,
      
      // Column 1: Panel E (Stacked Details)
      [
        #panel-label("E")
        #v(-10pt)
        #panel-title("Tracer concentration at right MCA")
        #v(-0pt)
        #grid(
          columns: (auto, 1fr), gutter: 2pt, align: horizon,
          
          y-label("high PVS \n flow ", nudge: 0.2cm),
          crop-image("../plots/modelA-strongVM/modelA-strongVM_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png", 0pt, 0pt, 17pt, 10pt),
          y-label("high PVS \n flow - dilated", nudge: -0.2cm),
          crop-image("../plots/LargePVS/LargePVS_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png", 10pt, 0pt, 7pt, 10pt),
                    v(0pt),
          x-label("concentration (mmol/l)", nudge: -1pt)
        ),     
      ],
      
      // Column 2: Panel F (Line Chart)
      [
        #panel-label("F")
        #v(-20pt)
        #panel-title("First-time arrival")
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
        #v(-25pt)
        #panel-title("Total tracer concentration")
        #v(-5pt)
        #cropped-image("../plots/comparisons/modelA_LargePVSA_modelA-strongVM_LargePVS/modelA_LargePVSA_modelA-strongVM_LargePVS_tot_bar.png", 5.1cm, 100%, 0cm, -0.1cm, width: 100%)
      ],
      
      // Column 2: Panel H (Line Charts Grid)
      [
        #panel-label("H")
        #v(-15pt)
        #panel-title("Mean tracer concentration")
        #v(-8pt)
        #cropped-image("../plots/comparisons/modelA_LargePVSA_modelA-strongVM_LargePVS/modelA_LargePVSA_modelA-strongVM_LargePVS_conc.png", 5.1cm, 100%, 0cm, 0.1cm, width: 110%)
      ],
    ),
    v(-30pt)
  )
)