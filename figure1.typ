#set page(width: 21cm, height: auto, margin: 0.5cm)
#set text(font: "Helvetica", size: 8pt)

#let panel-label(content) = align(top + left)[#text(weight: "bold", size: 12pt)[#content]]
#let panel-title(content) = align(center)[#text(weight: "bold", size: 8pt)[#content]]

// Main figure container
#figure(
  grid(
    columns: 1,
    row-gutter: 0pt, // Space between top, middle, and bottom rows

    // ==========================================
    // TOP ROW: Panels A, B, and C
    // ==========================================
    grid(
      columns: (2.5fr, 2.2fr, 1.5fr),
      gutter: 0pt,
      [
        #panel-label("A") 
        #v(-18pt)
        #image("paper/figures/Brain-PVS-callouts.png", width: 99%)
      ],
      [
        #panel-label("B")
        #v(-25pt) 
        #image("paper/figures/brain_mesh_vessels.png", width: 100%)
      ],
      [
        #panel-label("C") 
        #v(-15pt)
        #image("paper/figures/injection.png", width: 110%)
      ]
    ),

    // ==========================================
    // MIDDLE ROW: Panel D (Timeview 3D brains + colorbars)
    // ==========================================
    [
      #panel-label("D")
      #v(-40pt)
      #grid(
        columns: (1fr, 1fr, 1fr, 1fr, auto), // 4 brains, 1 column for legends
        align: horizon, // Vertically centers the brains and legends
        gutter: 0pt,
        image("plots/modelA/timeview3D/timeview3D_volume_3600.png", width: 100%),
        image("plots/modelA/timeview3D/timeview3D_volume_21600.png", width: 100%),
        image("plots/modelA/timeview3D/timeview3D_volume_43200.png", width: 100%),
        image("plots/modelA/timeview3D/timeview3D_volume_86400.png", width: 100%),
        
        // Nested grid specifically for the stacked colorbars and time labels
        grid(
          columns: (auto, auto),
          row-gutter: 0pt,
          column-gutter: 0pt,
          align: horizon,
          image("plots/modelA/timeview3D/colorbar_3600.png", height: 25pt),  [ 1 h ],
          image("plots/modelA/timeview3D/colorbar_21600.png", height: 25pt), [ 6 h ],
          image("plots/modelA/timeview3D/colorbar_43200.png", height: 25pt), [ 12 h ],
          image("plots/modelA/timeview3D/colorbar_86400.png", height: 25pt), [ 24 h ]
        )
      )
    ],

    // ==========================================
    // BOTTOM ROW: Panel E (Left) and Panels F & G (Right)
    // ==========================================
    grid(
      columns: (2.9fr, 1.2fr), // Left side is wide, right side is narrow
      gutter: 0pt,
      
      // Panel E
      [
        #panel-label("E") 
        #v(-13pt)
        #image("plots/modelA/modelA_overview_1-6-12-24.png", width: 100%)
      ],
      
      // Panels F & G stacked vertically
      grid(
        columns: 1,
        row-gutter: 0pt,
        [
          #panel-label("F")
          #v(-20pt)
          #panel-title("Mean tracer concentration")
          #v(-5pt)
          #image("plots/modelA/modelA_mean_conc.png", width: 100%)
        ],
        [
          #panel-label("G")
          #v(-20pt)
          #panel-title("Total tracer content")
          #v(-5pt)
          #image("plots/modelA/modelA_total_conc.png", width: 100%)
        ]
      )
    )
  )
)