#import "utils.typ": setup-figure, panel-label, panel-title, inline-header, cropped-image, y-label, x-label
#show: setup-figure

#let crop-image(path, t, r, b, l, bwidth: 100%, ..args) = {
  box(
    clip: true,
    width : bwidth,
    inset: (top: -t, right: -r, bottom: -b, left: -l),
    image(path, ..args)
  )
}


#let annotated-baseline(file) = {
  box(width: 100%, {
    // The image acts as the background layer
    crop-image(file, 0pt, 0pt, 20pt, 10pt)    

    // CSF Label & Line
    place(top + left, dx: 3%, dy: 25%)[#text(fill: white, size: 6pt)[CSF]]
    place(top + left, dx: 6%, dy: 35%)[#line(end: (3%, 8pt), stroke: white + 0.5pt)]
    
    // PVS Label & Line
    place(top + left, dx: 5%, dy: 83%)[#text(fill: white, size: 6pt)[PVS]]
    place(top + left, dx: 7%, dy: 81%)[#line(end: (5%, -10pt), stroke: white + 0.5pt)]
    
    // Artery Label & Line
    place(top + left, dx: 18%, dy: 30%)[#text(fill: white, size: 6pt)[artery]]
    place(top + left, dx: 21%, dy: 41%)[#line(end: (-21pt, 9pt), stroke: white + 0.5pt)]
    
    // Pia Label & Line
    place(top + left, dx: 12%, dy: 86%)[#text(fill: white, size: 6pt)[pia]]
    place(top + left, dx: 14%, dy: 85%)[#line(end: (7pt, -8pt), stroke: white + 0.5pt)]
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
        y-label("baseline", nudge: 0.2cm),
        cropped-image("../plots/modelA/modelA_conc_at_label_annotated.png", 3cm, 100%, 0cm, 0.15cm),
        
        // Label nudged by -0.55cm to match the image's move-y
        y-label("high PVS flow", nudge: -0.55cm),
        cropped-image("../plots/modelA-strongVM/modelA-strongVM_conc_at_label.png", 3.2cm, 100%, 0cm, -0.75cm)
      )
      #v(-10pt)
    ],

    // ------------------------------------------
    // MIDDLE BLOCK: 2x2 Grid (B & C top, E & D bottom)
    // ------------------------------------------
    grid(
      columns: (1fr, 1.1fr), 
      column-gutter: 1pt,
      row-gutter: 0pt,      
      
      // ROW 1, COLUMN 1: Panel B
      [
        #panel-label("B")
        #v(-15pt) #panel-title("Tracer first-time arrival") 
        #v(-5pt)
        #image("../plots/comparisons/modelA_modelA-strongVM_modelA-PVS-disp/modelA_modelA-strongVM_modelA-PVS-disp_fta.png", width: 100%)
      ],
      // ROW 1, COLUMN 2: Panel C
      [
        #panel-label("C")
        #v(-20pt) #panel-title("MCA-R") 
        #v(-10pt)
        #grid(
          columns: (auto, 1fr), gutter: 5pt, align: horizon,
          y-label("baseline", nudge: 0.2cm),
          annotated-baseline("../plots/modelA/modelA_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png"),
          
          y-label("high PVS flow", nudge: -0.20cm),
          crop-image("../plots/modelA-strongVM/modelA-strongVM_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png", 10pt, 0pt, 8pt, 10pt),
          v(0pt),
          x-label("concentration (mmol/l)", nudge: -3pt)
        )
      ],

      // ROW 2, COLUMN 1: Panel E
      [
        #panel-label("E")
        #v(-18pt) #panel-title("Total PVS tracer ")
        #v(-10pt)
        #grid(
          columns: (1fr, 1fr),
          gutter: 0pt,
          [
            #panel-title("baseline")
            #v(-8pt)
            #image("../plots/modelA/modelA_ridgeline_pvs_total_smoothed.png", width: 100%)
          ],
          [
            #panel-title("high PVS flow")
            #v(-8pt)
            #image("../plots/modelA-strongVM/modelA-strongVM_ridgeline_pvs_total_smoothed.png", width: 100%)
          ]
        )
      ],

      // ROW 2, COLUMN 2: Panel D
      [
        #panel-label("D")
        #v(-18pt) #panel-title("MCA2-R")
        #v(-10pt)
        #grid(
          columns: (auto, 1fr), gutter: 5pt, align: horizon,
          y-label("baseline", nudge: 0.2cm),      
          crop-image("../plots/modelA/modelA_2+4+6+8_MCA2-R_0.2+8.0+8.0+5.0_details.png", 0pt, 0pt, 20pt, 10pt),
          
          y-label("high PVS flow", nudge: -0.22cm), 
          crop-image("../plots/modelA-strongVM/modelA-strongVM_2+4+6+8_MCA2-R_0.2+8.0+8.0+5.0_details.png", 10pt, 0pt, 10pt, 10pt),
          v(0pt),
          x-label("concentration (mmol/l)", nudge: -3pt)
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
        #v(-15pt) #panel-title("Baseline") #v(-10pt)
        #cropped-image("../plots/modelA/modelA_overview_4-6.png", 9cm, 150%, 0cm, 0cm)
      ],
      [
        #panel-label("G")
        #v(-15pt) #panel-title("High PVS flow") #v(-10pt)
        #cropped-image("../plots/modelA-strongVM//modelA-strongVM_overview_4-6.png", 9cm, 150%, 0cm, 0cm)
      ],
      [
        #panel-label("H")
        #v(-15pt) #panel-title("High PVS dispersion") #v(-10pt)
        #cropped-image("../plots/modelA-PVS-disp/modelA-PVS-disp_overview_4-6.png", 9cm, 150%, 0cm, 0cm)
      ],
      [
        // Clever hack for the legend!
        #panel-label("")
        #v(-15pt) #panel-title("") #v(-10pt)
        #cropped-image("../plots/modelA-PVS-disp/modelA-PVS-disp_overview_4-6.png", 9cm, 300%, -5.5cm, 0cm)
      ]
    ),
    v(-20pt)
  )
)