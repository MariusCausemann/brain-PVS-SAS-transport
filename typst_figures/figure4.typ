#import "utils.typ": setup-figure, panel-label, panel-title, inline-header, cropped-image
#show: setup-figure

// ==========================================
// MAIN FIGURE LAYOUT
// ==========================================
#figure(
  grid(
    columns: 1,
    row-gutter: 20pt, 

    grid(
      columns: (4fr, 1fr), 
      gutter: 1pt,
      
      // Left Column (A)
      [
        #panel-label("A")
        #v(-25pt) 
        #image("../plots/meshplots/labeled_arteries.png", width: 100%)
      ],

      // Right Column (B)
      [
        #panel-label("B")
        #v(5pt) 
        #align(center)[
          #image("../paper/figures/peristaltic_flow_vertical.png", width: 70%)
        ]
      ]
    ),


    grid(
      columns: (1.3fr, 1fr, 1.3fr), 
      gutter: 15pt,
      
      // Column 1: Panel C (Cardiac 3D)
      [
        #panel-label("C")
        #v(-20pt) #panel-title("Cardiac-driven flow")
        #image("../plots/pvs_flow_peristaltic/cardiac_pvs_oscillation/vel3D_black.png", width: 145%)
      ],

      // Column 2: Panels D & F (Stacked Histograms)
      grid(
        columns: 1,
        row-gutter: 0pt,
        [
          #panel-label("D")
          #v(-20pt) #panel-title("Cardiac flow velocities")
          #v(-15pt)
          #image("../plots/pvs_flow_peristaltic/cardiac_pvs_oscillation/cardiac_pvs_oscillation_velocity_histo_cell.png", width: 120%)
        ],
        [
          #panel-label("F")
          #v(-20pt) #panel-title("Vasomotion flow velocities")
          #v(-15pt)
          #image("../plots/pvs_flow_peristaltic/vasomotion-strong/vasomotion-strong_velocity_histo_cell.png", width: 120%)
        ]
      ),

      // Column 3: Panel E (Vasomotion 3D)
      [
        #panel-label("E")
        #v(-20pt) #panel-title("Vasomotion-driven flow")
        #image("../plots/pvs_flow_peristaltic/vasomotion-strong/vel3D_black.png", width: 145%)
      ]
    )
  )
)