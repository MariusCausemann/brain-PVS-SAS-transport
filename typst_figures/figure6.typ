#import "utils.typ": setup-figure, panel-label, panel-title, inline-header, cropped-image, y-label, crop-image, x-label
#show: setup-figure


#let detail-block(title, img-low, img-high) = {
  grid(
    columns: 1, row-gutter: 1pt,
    panel-title(title, alignment: left, dx: -6pt),
    v(-5pt),
    grid(
      columns: (auto, 1fr), gutter: 2pt, align: horizon,
      y-label("low permeability", nudge: 0.1cm),
      crop-image(img-low, 0pt, 0pt, 18pt, 10pt),
      
      y-label("high permeability", nudge: -0.1cm),
      crop-image(img-high, 8pt, 0pt, 8pt, 10pt),
    ),
    x-label("concentration (mmol/l)")
  )
}

#let callout(start-x, start-y, end-x, end-y, off-x: 0pt, off-y: 0pt) = {
        // 1. Draw the box (centered around the start coordinate)
        place(top + left, dx: start-x, dy: start-y)[
          #move(dx: -6pt, dy: -6pt)[
            #rect(width: 12pt, height: 12pt, stroke: 1.5pt + black, fill: none)
          ]
        ]
        place(top + left, dx: start-x, dy: start-y)[
          #move(dx: off-x, dy: off-y)[
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
    box(width: 100%, [
      
      // THE GRID LAYOUT
      #grid(
        columns: (1fr, 1fr), 
        column-gutter: 15pt,
        row-gutter: 10pt,

        grid(
          columns: 1, row-gutter: 10pt,
          [
            #panel-label("A")
            #v(-5pt)
            #image("../paper/figures/pvs_permeability.png", width: 90%)
          ],
          [
            #panel-label("B")
            #v(-50pt)
            #cropped-image("../plots/meshplots/arteries.png", 9cm, 100%, 0cm, 1.3cm, width: 132%)
          ]
        ),
        
        [
          #panel-label("D", nudge-dx: -22pt)
          #v(-15pt) 
          #grid(
            columns: 1, row-gutter: 1pt,
            detail-block("ACA-A3",
              "../plots/modelA-strongVM/modelA-strongVM_2+4+6+8_ACA-A3_0.2+2.0+12.0+10.0_details.png",
              "../plots/modelB2-100/modelB2-100_2+4+6+8_ACA-A3_0.2+2.0+12.0+10.0_details.png"
            ),
            detail-block("ACA-A2",
              "../plots/modelA-strongVM/modelA-strongVM_2+4+6+8_ACA-A2_0.2+4.0+12.0+10.0_details.png",
              "../plots/modelB2-100/modelB2-100_2+4+6+8_ACA-A2_0.2+4.0+12.0+10.0_details.png"
            ),
            detail-block("MCA-R",
              "../plots/modelA-strongVM/modelA-strongVM_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png",
              "../plots/modelB2-100/modelB2-100_1+2+3+4_MCA-R_0.2+5.0+20.0+20.0_details.png"
            )
          )
          #v(-25pt),
        ],
        // --- ROW 2 ---
        [
          #panel-label("C", nudge-dy: -10pt)
          #v(-15pt)
          #detail-block("BA",
            "../plots/modelA-strongVM/modelA-strongVM_1+2+3+4_BA_3.0+40.0+40.0+10.0_details.png",
            "../plots/modelB2-100/modelB2-100_1+2+3+4_BA_3.0+40.0+40.0+10.0_details.png"
          )
        ],
        
        [
          #panel-label("")
          #v(-15pt) 
          #detail-block("MCA-L",
            "../plots/modelA-strongVM/modelA-strongVM_1+2+3+4_MCA-L_0.2+5.0+20.0+20.0_details.png",
            "../plots/modelB2-100/modelB2-100_1+2+3+4_MCA-L_0.2+5.0+20.0+20.0_details.png"
          )
        ]
      )
      
      #callout(29.5%, 51%, 4cm, -6cm,   off-x: 6pt, off-y: -6pt) 
      #callout(27%, 54%, 4.3cm, -5cm, off-x: 6pt, off-y: 0pt) 
      #callout(21%, 60%, 5.3cm, -1.6cm,   off-x: 6pt, off-y: 6pt) 
      #callout(33.4%, 60%, 3.1cm, 2.8cm,   off-x: 6pt, off-y: 6pt) 
      #callout(25%, 66.5%, -4.5cm, 1.7cm,  off-x: -6pt, off-y: 6pt)
    ]),

    // ------------------------------------------
    // BOTTOM SECTION: Panels E and F
    // ------------------------------------------
    [
      #panel-label("E")
      #v(-16pt)
      #image("../plots/comparisons/modelA_modelA-strongVM_modelB2-10_modelB2-100/modelA_modelA-strongVM_modelB2-10_modelB2-100_conc.png", width: 100%)
    ],
    [
      #panel-label("F")
      #v(-16pt)
      #image("../plots/comparisons/modelA-strongVM_modelB2-10_modelB2-100/modelA-strongVM_modelB2-10_modelB2-100.png", width: 100%)
    ],
    v(-20pt)
  )
)