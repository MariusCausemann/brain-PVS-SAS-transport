#set page(width: 21cm, height: auto, margin: 1cm)
#set text(font: "Nimbus Sans", size: 8pt)

// ==========================================
// HELPER FUNCTIONS
// ==========================================
#let panel-label(content) = align(top + left)[#text(weight: "bold", size: 12pt)[#content]]
#let panel-title(content) = align(center)[#text(weight: "bold", size: 8pt)[#content]]

#let y-label(content, nudge: 0cm) = align(center + horizon)[
  #move(dy: nudge)[
    #rotate(-90deg, reflow: true)[#text(size: 6.5pt)[#content]]
  ]
]

#let cropped-image(file, box-height, zoom, move-x, move-y) = {
  box(
    width: 100%, 
    height: box-height, 
    clip: true,
    place(
      dx: move-x, 
      dy: move-y, 
      image(file, width: zoom)
    )
  )
}

#let annotated-baseline(file) = {
  box(width: 100%, {
    // The image acts as the background layer
    cropped-image(file, 2.3cm, 100%, 0cm, 0.49cm)    
    
    // CSF Label & Line
    place(top + left, dx: 5%, dy: 40%)[#text(fill: white, size: 6pt)[CSF]]
    place(top + left, dx: 7%, dy: 48%)[#line(end: (5%, 7pt), stroke: white + 0.5pt)]
    
    // PVS Label & Line
    place(top + left, dx: 5%, dy: 83%)[#text(fill: white, size: 6pt)[PVS]]
    place(top + left, dx: 7%, dy: 78%)[#line(end: (7.5%, -8pt), stroke: white + 0.5pt)]
    
    // Artery Label & Line
    place(top + left, dx: 21%, dy: 40%)[#text(fill: white, size: 6pt)[artery]]
    place(top + left, dx: 24%, dy: 47%)[#line(end: (-20pt, 12pt), stroke: white + 0.5pt)]
    
    // Pia Label & Line
    place(top + left, dx: 17%, dy: 86%)[#text(fill: white, size: 6pt)[pia]]
    place(top + left, dx: 19%, dy: 85%)[#line(end: (4pt, -5pt), stroke: white + 0.5pt)]
  })
}

// ==========================================
// MAIN FIGURE LAYOUT
// ==========================================
#figure(
  grid(
    columns: 1,
    row-gutter: 0pt, 

    // ------------------------------------------
    // TOP BLOCK: Panel A 
    // ------------------------------------------
    [
      #panel-label("A")
      #v(-18pt)
      #grid(
        columns: (auto, 1fr), 
        row-gutter: 5pt,
        align: horizon, 
        // Label nudged by 0.15cm to match the image's move-y
        y-label("baseline", nudge: 0.15cm),
        cropped-image("plots/modelA/modelA_conc_at_label_annotated.png", 3cm, 100%, 0cm, 0.15cm),
        
        // Label nudged by -0.55cm to match the image's move-y
        y-label("high PVS flow", nudge: -0.55cm),
        cropped-image("plots/modelA-strongVM/modelA-strongVM_conc_at_label.png", 3cm, 100%, 0cm, -0.55cm)
      )
      #v(-12pt)
    ],

    // ------------------------------------------
    // MIDDLE BLOCK: 2x2 Grid (B & C top, E & D bottom)
    // ------------------------------------------
    grid(
      columns: (1fr, 1fr), 
      column-gutter: 0pt,
      row-gutter: 0pt,      
      
      // ROW 1, COLUMN 1: Panel B
      [
        #panel-label("B")
        #v(-12pt)
        #image("plots/comparisons/modelA_modelA-strongVM_modelA-PVS-disp/modelA_modelA-strongVM_modelA-PVS-disp_fta.png", width: 100%)
      ],

      // ROW 1, COLUMN 2: Panel C
      [
        #panel-label("C")
        #v(-20pt) #panel-title("MCA-R") 
        #v(-18pt)
        #grid(
          columns: (auto, 1fr), gutter: 5pt, align: horizon,
          y-label("baseline", nudge: 0.49cm),
          annotated-baseline("plots/modelA/modelA_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png"),
          
          y-label("high PVS flow", nudge: -0.22cm),
          cropped-image("plots/modelA-strongVM/modelA-strongVM_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png", 2.3cm, 100%, 0cm, -0.22cm)
        )
      ],

      // ROW 2, COLUMN 1: Panel E
      [
        #panel-label("E")
        #v(-20pt) #panel-title("total PVS tracer content")
        #v(-10pt)
        #grid(
          columns: (1fr, 1fr),
          gutter: 5pt,
          [
            #panel-title("baseline")
            #image("plots/modelA/modelA_ridgeline_pvs_total_smoothed.png", width: 100%)
          ],
          [
            #panel-title("high PVS flow")
            #image("plots/modelA-strongVM/modelA-strongVM_ridgeline_pvs_total_smoothed.png", width: 100%)
          ]
        )
      ],

      // ROW 2, COLUMN 2: Panel D
      [
        #panel-label("D")
        #v(-10pt) #panel-title("MCA2-R")
        #v(-18pt)
        #grid(
          columns: (auto, 1fr), gutter: 5pt, align: horizon,
          y-label("baseline", nudge: 0.49cm),      
          cropped-image("plots/modelA/modelA_2+4+6+8_MCA2-R_0.2+8.0+8.0+5.0_details.png", 2.3cm, 100%, 0cm, 0.49cm),
          
          y-label("high PVS flow", nudge: -0.22cm), 
          cropped-image("plots/modelA-strongVM/modelA-strongVM_2+4+6+8_MCA2-R_0.2+8.0+8.0+5.0_details.png", 2.3cm, 100%, 0cm, -0.22cm)
        )
      ]
    ),

    // ------------------------------------------
    // BOTTOM BLOCK: Panels F, G, H
    // ------------------------------------------
    grid(
      columns: (1fr, 1fr, 1fr, 0.5fr), 
      gutter: 0pt, 
      [
        #panel-label("F")
        #v(-20pt) #panel-title("baseline") #v(-2pt)
        #cropped-image("plots/modelA/modelA_overview_4-6.png", 9cm, 150%, 0cm, 0cm)
      ],
      [
        #panel-label("G")
        #v(-20pt) #panel-title("high PVS flow") #v(-2pt)
        #cropped-image("plots/modelA-strongVM//modelA-strongVM_overview_4-6.png", 9cm, 150%, 0cm, 0cm)
      ],
      [
        #panel-label("H")
        #v(-20pt) #panel-title("high PVS dispersion") #v(-2pt)
        #cropped-image("plots/modelA-PVS-disp/modelA-PVS-disp_overview_4-6.png", 9cm, 150%, 0cm, 0cm)
      ],
      [
        // Clever hack for the legend!
        #panel-label("")
        #v(-20pt) #panel-title("") #v(-2pt)
        #cropped-image("plots/modelA-PVS-disp/modelA-PVS-disp_overview_4-6.png", 9cm, 300%, -5.5cm, 0cm)
      ]
    )
  )
)