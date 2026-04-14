#import "utils.typ": setup-figure, panel-label, panel-title
#show: setup-figure

// Main figure container (Part 2)
#figure(
  grid(
    columns: 1,
    row-gutter: 10pt, 

    // ==========================================
    // TOP ROW: Panels C and D
    // ==========================================
    grid(
      columns: (0.8fr, 3fr), // D gets much more space to fit the 4 brains + legend
      gutter: 5pt,
      [
        #panel-label("A") 
        #v(-15pt)
        #image("../paper/figures/injection.png", width: 105%)
        #v(-20pt)
      ],
      [
        #panel-label("B")
        #v(-10pt)
        #grid(
          columns: (1fr, 1fr, 1fr, 1fr, auto), // 4 brains, 1 column for legends
          align: horizon, 
          gutter: 0pt,
          image("../plots/modelA/timeview3D/timeview3D_volume_3600.png", width: 110%),
          image("../plots/modelA/timeview3D/timeview3D_volume_21600.png", width: 110%),
          image("../plots/modelA/timeview3D/timeview3D_volume_43200.png", width: 110%),
          image("../plots/modelA/timeview3D/timeview3D_volume_86400.png", width: 110%),
          
          // Nested grid for the stacked colorbars and time labels
          grid(
            columns: (auto, auto),
            row-gutter: 0pt,
            column-gutter: 0pt,
            align: horizon,
            image("../plots/modelA/timeview3D/colorbar_3600.png", height: 30pt),  [ 1 h ],
            image("../plots/modelA/timeview3D/colorbar_21600.png", height: 30pt), [ 6 h ],
            image("../plots/modelA/timeview3D/colorbar_43200.png", height: 30pt), [ 12 h ],
            image("../plots/modelA/timeview3D/colorbar_86400.png", height: 30pt), [ 24 h ]
          )
        )
      #v(-20pt)
      ]
    ),

    // ==========================================
    // BOTTOM ROW: Panel E (Left) and Panels F & G (Right)
    // ==========================================
    grid(
      columns: (2.9fr, 1.2fr), 
      gutter: 1pt, // Added a small gutter so E doesn't touch F/G
      
      // Panel E
      [
        #panel-label("C") 
        #v(-13pt)
        #image("../plots/modelA/modelA_overview_1-6-12-24.png", width: 100%)
      ],
      
      // Panels F & G stacked vertically
      grid(
        columns: 1,
        row-gutter: 0pt,
        [
          #panel-label("D")
          #v(-20pt)
          #panel-title("Mean tracer concentration")
          #v(-8pt)
          #image("../plots/modelA/modelA_mean_conc.png", width: 100%)
        ],
        [
          #panel-label("E")
          #v(-20pt)
          #panel-title("Total tracer content")
          #v(-8pt)
          #image("../plots/modelA/modelA_total_conc.png", width: 100%)
        ]
      )
    ),
    v(-20pt)
  )
)