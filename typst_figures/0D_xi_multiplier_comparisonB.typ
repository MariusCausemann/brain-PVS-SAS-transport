#import "utils.typ": setup-figure, panel-label, panel-title, cropped-image, crop-image
#show: setup-figure

#let titlefontsize = 10pt

#figure(
  grid(
    columns: 1,
    row-gutter: 15pt, 
    v(-40pt),
    grid(
      columns: (1.2fr, 1.2fr, 0.4fr), 
      gutter: 0pt, 
      [
        #panel-label("A")
        #v(-15pt) #panel-title("Overview (with 0D)", fontsize: titlefontsize) 
        #v(-10pt)
        #cropped-image("../plots/modelA-strongVM/modelA-strongVM_overview_9-12-24.png", 11cm, 130%, 0cm, 0cm)
      ],
      [
        #panel-label("B")
        #v(-15pt) #panel-title("Overview (without 0D)", fontsize: titlefontsize) 
        #v(-10pt)
        #cropped-image("../plots/modelA-strongVM-root100/modelA-strongVM-root100_overview_9-12-24.png", 11cm, 130%, 0cm, 0cm)
      ],
      [
        #panel-label("")
        #v(-15pt) #panel-title("") #v(-10pt)
        #cropped-image("../plots/modelA-strongVM/modelA-strongVM_overview_9-12-24.png", 11cm, 390%, -8.2cm, 0cm)
      ]
    ),
    v(-100pt),
    grid(
      columns: (1fr), 
      gutter: 0pt, 
      [
        #panel-label("C")
        #v(-15pt) #panel-title("Concentration difference (with 0D - without 0D)", fontsize: titlefontsize) 
        #v(-5pt)
        #image("../plots/comparisons/modelA-strongVM_vs_modelA-strongVM-root100_diff_1-6-9-12-24.png", height: 9.0cm)
      ],

    )
  ),
)
#v(-30pt)
