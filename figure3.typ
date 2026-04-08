#set page(width: 21cm, height: auto, margin: 1cm)
#set text(font: "Nimbus Sans", size: 8pt)


// ==========================================
// HELPER FUNCTIONS
// ==========================================
#let panel-label(content) = align(top + left)[#text(weight: "bold", size: 12pt)[#content]]
#let panel-title(content) = align(center)[#text(weight: "bold", size: 8pt)[#content]]

// ==========================================
// MAIN FIGURE LAYOUT
// ==========================================
#figure(
  grid(
    columns: 1,
    row-gutter: 20pt, 

    // ------------------------------------------
    // TOP BLOCK: Panel A (Left) | Panel B (Right)
    // ------------------------------------------
    grid(
      columns: (1.6fr, 1fr), // Proportions based on your 524px vs 320px widths
      gutter: 15pt,
      
      // Left Column (A)
      [
        #panel-label("A")
        #image("plots/meshplots/labeled_arteries.png", width: 100%)
      ],

      // Right Column (B - Stacked 3D plot and histogram)
      [
        #panel-label("B")
        #v(-5pt) #panel-title("pressure-driven flow")
        #grid(
          columns: 1, 
          row-gutter: 10pt,
          image("plots/pvs_flow_prod/sas_flow-arteries-2/vel3D.png", width: 100%),
          image("plots/pvs_flow_prod/sas_flow-arteries/sas_flow-arteries_velocity_histo_cell.png", width: 100%)
        )
      ]
    ),

    // ------------------------------------------
    // BOTTOM BLOCK: Panels C | D (Stacked) | E | F
    // ------------------------------------------
    grid(
      columns: (0.3fr, 1fr, 1.4fr, 1.4fr), // C is narrow, E/F are wide
      gutter: 15pt,
      
      // Column 1 (C - Narrow Vertical Diagram)
      [
        #panel-label("C")
        // Added top margin to push it down slightly to align with titles in D/E/F
        #v(15pt) 
        #image("paper/figures/peristaltic_flow_vertical.png", width: 100%)
      ],

      // Column 2 (D - Stacked Histograms)
      [
        #panel-label("D")
        #grid(
          columns: 1,
          row-gutter: 10pt,
          [
            #v(-5pt) #panel-title("cardiac-driven flow")
            #image("plots/pvs_flow_peristaltic/cardiac_pvs_oscillation/cardiac_pvs_oscillation_velocity_histo_cell.png", width: 100%)
          ],
          [
            #v(5pt) // Extra spacing before the second title
            #panel-title("vasomotion-driven flow")
            #image("plots/pvs_flow_peristaltic/vasomotion-strong/vasomotion-strong_velocity_histo_cell.png", width: 100%)
          ]
        )
      ],

      // Column 3 (E - Cardiac 3D)
      [
        #panel-label("E")
        #v(-5pt) #panel-title("cardiac-driven flow")
        #image("plots/pvs_flow_peristaltic/cardiac_pvs_oscillation/vel3D.png", width: 100%)
      ],

      // Column 4 (F - Vasomotion 3D)
      [
        #panel-label("F")
        #v(-5pt) #panel-title("vasomotion-driven flow")
        #image("plots/pvs_flow_peristaltic/vasomotion-strong/vel3D.png", width: 100%)
      ]
    )
  )
)