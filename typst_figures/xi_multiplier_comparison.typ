#import "utils.typ": setup-figure, panel-label, panel-title
#show: setup-figure

#let titlefontsize = 10pt
#let subtitlefontsize = 8pt

#figure(
  grid(
    columns: 1,
    gutter: 15pt, 

    // ---------------------------------------------------------
    // ROW 1: MACRO-ROBUSTNESS (Total Content)
    // ---------------------------------------------------------
    [
      #panel-label("A")
      #v(-15pt) 
      #panel-title("Total tracer content", fontsize: titlefontsize)
      #v(-5pt)
      #grid(
        columns: (1fr, 1fr, 1fr, 1fr),
        gutter: 0pt,
        [#image("../plots/modelA-strongVM-root1/modelA-strongVM-root1_total_conc.png", width: 105%) #v(-5pt) #align(center)[#text(size: subtitlefontsize)[1x multiplier]]],
        [#image("../plots/modelA-strongVM-root10/modelA-strongVM-root10_total_conc.png", width: 105%) #v(-5pt) #align(center)[#text(size: subtitlefontsize)[10x multiplier]]],
        [#image("../plots/modelA-strongVM-root100/modelA-strongVM-root100_total_conc.png", width: 105%) #v(-5pt) #align(center)[#text(size: subtitlefontsize)[100x multiplier]]],
        [#image("../plots/modelA-strongVM-root1000/modelA-strongVM-root1000_total_conc.png", width: 105%) #v(-5pt) #align(center)[#text(size: subtitlefontsize)[1000x multiplier]]]
      )
    ],

    // ---------------------------------------------------------
    // ROW 2: KINETIC ROBUSTNESS (First-time Arrival)
    // ---------------------------------------------------------
    grid(
      columns: 1,
      gutter: 0pt,
      [
        #panel-label("B")
        #v(-15pt) 
        #panel-title("First-time arrival", fontsize: titlefontsize)
        #v(-8pt)
        #align(center)[
        #image("../plots/comparisons/modelA-strongVM-root1_modelA-strongVM-root10_modelA-strongVM-root100_modelA-strongVM-root1000/modelA-strongVM-root1_modelA-strongVM-root10_modelA-strongVM-root100_modelA-strongVM-root1000_fta.png", width: 80%)
        ]
      ]
    ),

    // ---------------------------------------------------------
    // ROW 3: LOCAL ARTIFACT (Leaf Concentrations)
    // ---------------------------------------------------------
    [
      #panel-label("C")
      #v(-15pt) 
      #panel-title("Arterial leaf node concentrations", fontsize: titlefontsize)
      #v(-5pt)
      #grid(
        columns: (1fr, 1fr, 1fr, 1fr),
        gutter: 0pt,
        [#image("../plots/modelA-strongVM-root1/modelA-strongVM-root1_leaf_conc.png", width: 110%) #v(-5pt) #align(center)[#text(size: subtitlefontsize)[1x multiplier]]],
        [#image("../plots/modelA-strongVM-root10/modelA-strongVM-root10_leaf_conc.png", width: 110%) #v(-5pt) #align(center)[#text(size: subtitlefontsize)[10x multiplier]]],
        [#image("../plots/modelA-strongVM-root100/modelA-strongVM-root100_leaf_conc.png", width: 110%) #v(-5pt) #align(center)[#text(size: subtitlefontsize)[100x multiplier]]],
        [#image("../plots/modelA-strongVM-root1000/modelA-strongVM-root1000_leaf_conc.png", width: 110%) #v(-5pt) #align(center)[#text(size: subtitlefontsize)[1000x multiplier]]]
      )
      #v(-20pt)
    ]
  )
)