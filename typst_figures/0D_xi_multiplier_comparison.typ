#import "utils.typ": setup-figure, panel-label, panel-title, cropped-image, crop-image
#show: setup-figure

#let titlefontsize = 10pt

#figure(
  grid(
    columns: 1,
    row-gutter: 15pt, 

// ---------------------------------------------------------
    // ROW 1: NEW 0D MODEL (A | B | C)
    // ---------------------------------------------------------
    grid(
      columns: (1fr, 1fr, 1fr),
      column-gutter: 10pt,
      
      [
        #panel-label("A") #v(-15pt) 
        #panel-title("Tracer content (0D model)", fontsize: titlefontsize)
        #v(-10pt)
        #image("../plots/modelA-strongVM/modelA-strongVM_total_conc.png", width: 100%)
      ],
      [
        #panel-label("B") #v(-15pt) 
        #panel-title("Tracer content (without 0D model)", fontsize: titlefontsize)
        #v(-10pt)
        #image("../plots/modelA-strongVM-root100/modelA-strongVM-root100_total_conc.png", width: 100%)
      ],
      [
        #panel-label("C") #v(-15pt) 
        #panel-title("0D concentrations (0D model)", fontsize: titlefontsize)
        #crop-image("../plots/modelA-strongVM/modelA-strongVM_tip_conc.png",0pt,0pt,13pt,0pt)

      ]
    ),

// ---------------------------------------------------------
    // ROW 2: OLD MODEL COMPARISON (E 2/3 width | D & H stacked 1/3 width)
    // ---------------------------------------------------------
    grid(
      columns: (2fr, 1fr),
      column-gutter: 10pt,
      
      // Panel E: First-time arrival (Taking 2/3 width)
      [
        #panel-label("F") #v(-20pt) 
        #panel-title("First-time tracer arrival", fontsize: titlefontsize)
        #v(-5pt)
        #image("../plots/comparisons/modelA-strongVM_modelA-strongVM-root100/modelA-strongVM_modelA-strongVM-root100_fta.png", width: 100%)
      ],

      // Right Column: Old model metrics (Stacked to fill the 1/3 width)
      grid(
        columns: 1,
        row-gutter: 0pt,
        [
          #v(-35pt)
          #panel-label("D")  
          #panel-title("Leaf concentrations (0D model)", fontsize: titlefontsize)
          #v(-12pt)
          #crop-image("../plots/modelA-strongVM/modelA-strongVM_leaf_conc.png",0pt,0pt,13pt,0pt)
        ],
        [
          #v(-15pt)
          #panel-label("E")
          #panel-title("Leaf concentrations (without 0D model)", fontsize: titlefontsize)
          #v(-12pt) 
          #image("../plots/modelA-strongVM-root100/modelA-strongVM-root100_leaf_conc.png", width: 100%)

        ]
      )
    ),

    // ---------------------------------------------------------
    // BOTTOM SECTION: Spatial Overview (Side-by-Side)
    // ---------------------------------------------------------
    v(-40pt),
    grid(
      columns: (1.2fr, 1.2fr, 0.4fr), 
      gutter: 0pt, 
      [
        #panel-label("G")
        #v(-15pt) #panel-title("Overview (with 0D)", fontsize: titlefontsize) 
        #v(-10pt)
        #cropped-image("../plots/modelA-strongVM/modelA-strongVM_overview_9-12-24.png", 11cm, 130%, 0cm, 0cm)
      ],
      [
        #panel-label("H")
        #v(-15pt) #panel-title("Overview (without 0D)", fontsize: titlefontsize) 
        #v(-10pt)
        #cropped-image("../plots/modelA-strongVM-root100/modelA-strongVM-root100_overview_9-12-24.png", 11cm, 130%, 0cm, 0cm)
      ],
      [
        #panel-label("")
        #v(-15pt) #panel-title("") #v(-10pt)
        #cropped-image("../plots/modelA-strongVM/modelA-strongVM_overview_9-12-24.png", 11cm, 390%, -8.2cm, 0cm)
      ]
    )
  ),
)
#v(-90pt)
