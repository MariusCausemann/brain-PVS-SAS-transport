#import "utils.typ": setup-figure, panel-label, panel-title
#show: setup-figure

#let titlefontsize = 10pt

#figure(
  grid(
    columns: 1,
    row-gutter: 10pt, 

    grid(
      columns: (2.4fr, 1fr), 
      gutter: 10pt,
      
      // Left Column: Panel A
      [
        #panel-label("A")
        #v(-20pt) 
        #panel-title("Physical mechanisms of intracranial solute transport", fontsize: titlefontsize)
        #image("../paper/figures/transport_mechanisms.png", width: 98%)
      ],
      
      // Right Column: Panels B & C stacked
      grid(
        columns: 1,
        row-gutter: 0pt,
        [
          #panel-label("B")
          #v(-20pt) 
          #panel-title("Intracranial 3D domains", fontsize: titlefontsize)
          #v(-15pt) 
          #image("../plots/meshplots/subdomains.png", width: 115%)
        ],
        [
          #panel-label("C")
          #v(-10pt) 
          #panel-title("Computational boundaries",fontsize: titlefontsize)
          #v(-25pt) 
          #image("../plots/meshplots/boundaries.png", width: 115%)
        ]
      )
    ),

    // BOTTOM SECTION: 1D Embedded Networks
    grid(
      columns: (1fr, 1fr, 1fr),
      gutter: 0pt,
      
      [
        #panel-label("D")
        #v(-15pt) 
        #panel-title("3D domains with PVS networks", fontsize: titlefontsize)
        #v(-15pt) 
        #image("../paper/figures/brain_mesh_vessels.png", width: 100%)
      ],
      [
        #panel-label("E")
        #v(-15pt) 
        #panel-title("Embedded arterial network", fontsize: titlefontsize)
        #v(-15pt) 
        #image("../paper/figures/arteries_red.png", width: 120%)
      ],
      [
        #panel-label("F")
        #v(-15pt) 
        #panel-title("Embedded venous network", fontsize: titlefontsize)
        #v(-15pt) 
        #image("../paper/figures/veines.png", width: 120%)
      ]
    ),
    v(-80pt)
  )
)