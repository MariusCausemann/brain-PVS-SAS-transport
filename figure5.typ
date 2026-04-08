#set page(width: 21cm, height: auto, margin: 1cm)
#set text(font: "Nimbus Sans", size: 8pt)

// ==========================================
// HELPER FUNCTIONS
// ==========================================

#let panel-label(content, nudge-dx: 0pt, nudge-dy: 0pt) = align(top + left)[
  #move(dx: nudge-dx, dy: nudge-dy)[#text(weight: "bold", size: 12pt)[#content]]
]

#let panel-title(content) = align(left)[#text(weight: "bold", size: 8pt)[#content]]

#let y-label(content, nudge: 0cm) = align(center + horizon)[
  #move(dy: nudge)[
    #rotate(-90deg, reflow: true)[#text(size: 6.5pt)[#content]]
  ]
]

#let cropped-image(file, box-height, zoom, move-x, move-y, width: 100%) = {
  box(
    width: width, 
    height: box-height, 
    clip: true,
    place(
      dx: move-x, 
      dy: move-y, 
      image(file, width: zoom)
    )
  )
}

#let detail-block(title, img-low, img-high) = {
  grid(
    columns: 1, row-gutter: 1pt,
    panel-title(title),
    v(-20pt),
    grid(
      columns: (auto, 1fr), gutter: 2pt, align: horizon,
      y-label("low permeability", nudge: 0.49cm),
      cropped-image(img-low, 2.3cm, 100%, 0cm, 0.49cm),
      
      y-label("high permeability", nudge: -0.22cm),
      cropped-image(img-high, 2.3cm, 100%, 0cm, -0.22cm)
    )
  )
}

#let callout(start-x, start-y, end-x, end-y, off-x: 0pt, off-y: 0pt) = {
        // 1. Draw the box (centered around the start coordinate)
        place(top + left, dx: start-x, dy: start-y)[
          #move(dx: -6pt, dy: -6pt)[
            #rect(width: 12pt, height: 12pt, stroke: 1.5pt + black, fill: none)
          ]
        ]
        // 2. Draw the line, nudged to start at the edge of the box
        place(top + left, dx: start-x, dy: start-y)[
          #move(dx: off-x, dy: off-y)[
            // We subtract the offset from the end destination so it still hits the exact same target!
            #line(end: (end-x - off-x, end-y - off-y), stroke: (paint: black, thickness: 1.5pt, dash: "dashed"))
          ]
        ]
      }

// ==========================================
// MAIN FIGURE LAYOUT
// ==========================================
#figure(
  grid(
    columns: 1,
    row-gutter: 0pt, 

    // ------------------------------------------
    // TOP SECTION: Canvas for grids and lines
    // ------------------------------------------
    box(width: 100%, [
      
      // THE GRID LAYOUT
      #grid(
        columns: (1fr, 1fr), 
        column-gutter: 15pt,
        row-gutter: 10pt,

        // --- ROW 1 ---
        grid(
          columns: 1, row-gutter: 10pt,
          [
            #panel-label("A")
            #v(-5pt)
            #image("paper/figures/pvs_permeability.png", width: 90%)
          ],
          [
            #panel-label("B")
            #v(-50pt)
            #cropped-image("plots/meshplots/arteries.png", 9cm, 100%, 0cm, 1.3cm, width: 132%)
          ]
        ),
        
        [
          #panel-label("D", nudge-dx: -12pt)
          #v(-15pt) 
          #grid(
            columns: 1, row-gutter: 0pt,
            detail-block("ACA-A3",
              "plots/modelA-strongVM/modelA-strongVM_2+4+6+8_ACA-A3_0.2+2.0+12.0+10.0_details.png",
              "plots/modelB2-100/modelB2-100_2+4+6+8_ACA-A3_0.2+2.0+12.0+10.0_details.png"
            ),
            detail-block("ACA-A2",
              "plots/modelA-strongVM/modelA-strongVM_2+4+6+8_ACA-A2_0.2+4.0+12.0+10.0_details.png",
              "plots/modelB2-100/modelB2-100_2+4+6+8_ACA-A2_0.2+4.0+12.0+10.0_details.png"
            ),
            detail-block("MCA-R",
              "plots/modelA-strongVM/modelA-strongVM_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png",
              "plots/modelB2-100/modelB2-100_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png"
            )
          )
          #v(-30pt),
        ],
        // --- ROW 2 ---
        [
          #panel-label("C", nudge-dy: -12pt)
          #v(-15pt)
          #detail-block("BA",
            "plots/modelA-strongVM/modelA-strongVM_1+2+3+4_BA_3.0+40.0+40.0+10.0_details.png",
            "plots/modelB2-100/modelB2-100_1+2+3+4_BA_3.0+40.0+40.0+10.0_details.png"
          )
        ],
        
        [
          #panel-label("")
          #v(-15pt) 
          #detail-block("MCA-L",
            "plots/modelA-strongVM/modelA-strongVM_1+2+3+4_MCA-L_0.2+5.0+20.0+20.0_details.png",
            "plots/modelB2-100/modelB2-100_1+2+3+4_MCA-L_0.2+5.0+20.0+20.0_details.png"
          )
        ]
      )
      
// We use off-x and off-y to push the start of the line 6pt toward its destination
      #callout(29.5%, 53%, 4cm, -6cm,   off-x: 6pt, off-y: -6pt) // Top: ACA-A3 (Line goes Right & Up)
      #callout(27%, 57%, 4.5cm, -5cm, off-x: 6pt, off-y: 0pt)  // Middle: ACA-A2 (Line goes Right)
      #callout(21%, 63%, 5.7cm, -2cm,   off-x: 6pt, off-y: 6pt)  // Bot-Right: MCA-L (Line goes Right & Down)
      #callout(33.4%, 63%, 3.3cm, 2.1cm,   off-x: 6pt, off-y: 6pt)  // Bot-Mid: MCA-R (Line goes Right & Down)
      #callout(25%, 70%, -4.3cm, 1cm,  off-x: -6pt, off-y: 6pt) // Bot-Left: BA (Line goes Left & Down)

    ]),

    // ------------------------------------------
    // BOTTOM SECTION: Panels E and F
    // ------------------------------------------
    [
      #panel-label("E")
      #v(-20pt)
      #image("plots/comparisons/modelA_modelA-strongVM_modelB2-10_modelB2-100/modelA_modelA-strongVM_modelB2-10_modelB2-100_conc.png", width: 100%)
    ],
    [
      #panel-label("F")
      #v(-20pt)
      #image("plots/comparisons/modelA-strongVM_modelB2-10_modelB2-100/modelA-strongVM_modelB2-10_modelB2-100.png", width: 100%)
    ]
  )
)