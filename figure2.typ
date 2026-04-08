#set page(width: 21cm, height: auto, margin: 1cm)
#set text(font: "Helvetica", size: 8pt)

// ==========================================
// HELPER FUNCTIONS
// ==========================================
#let panel-label(content) = align(top + left)[#text(weight: "bold", size: 12pt)[#content]]
#let panel-title(content) = align(center)[#text(weight: "bold", size: 8pt)[#content]]

// The Typst clipping hack
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

// ==========================================
// MAIN FIGURE LAYOUT
// ==========================================
#figure(
  grid(
    columns: 1,
    row-gutter: 0pt, 

    // ------------------------------------------
    // TOP BLOCK: (A, B) | C | (D, E)
    // ------------------------------------------
    grid(
      columns: (1fr, 1.85fr, 1fr), 
      gutter: 10pt,
      
      // Left Column (A & B)
      grid(
        columns: 1, row-gutter: 10pt,
        [
          #panel-label("A")
          #v(-15pt) 
          #panel-title("3D model of SAS and ventricles")
          #v(-15pt) 
          #image("plots/meshplots/subdomains.png", width: 100%)
        ],
        [
          #v(-25pt) 
          #panel-label("B")
          #v(-15pt) 
          #panel-title("boundaries of the CSF space")
          #v(-15pt) 
          #image("plots/meshplots/boundaries.png", width: 100%)
        ]
      ),

      // Middle Column (C)
      [
        #panel-label("C")
        #v(-15pt) 
        #panel-title("model illustration")
        #image("paper/figures/flow_and_dispersion.png", width: 110%)
      ],

      // Right Column (D & E)
      grid(
        columns: 1, row-gutter: 10pt,
        [
          #panel-label("D")
          #v(-15pt) 
          #panel-title("production-driven CSF flow")
          #v(-15pt) 
          #image("results/csf_flow/sas_flow/csf_v.png", width: 100%)
        ],
        [
          #panel-label("E")
          #v(-5pt) #panel-title("production-driven CSF velocities")
          #image("results/csf_flow/sas_flow/prod_velocity_histo.png", width: 100%)
        ]
      )
    ),

    // ------------------------------------------
    // MIDDLE BLOCK: Flow plots stacked over Histograms
    // ------------------------------------------
    grid(
      columns: (1fr, 1fr, 1fr, 1fr),
      gutter: 0pt,
      
      grid(
        columns: 1, row-gutter: 5pt,
        [
          #panel-label("F")
          #v(-5pt) 
          #panel-title("cardiac-driven CSF flow")
          #v(-20pt) 
          #image("results/csf_flow/cardiac_sas_flow/csf_v.png", width: 100%)
        ],
        [
          #panel-label("J")
          #v(-15pt) 
          #panel-title("cardiac flow velocities")
          #v(-5pt) 
          #image("results/csf_flow/cardiac_sas_flow/cardiac_velocity_histo.png", width: 100%)
        ]
      ),

      grid(
        columns: 1, row-gutter: 5pt,
        [
          #panel-label("G")
          #v(-5pt) 
          #panel-title("cardiac dispersion enhancement")
          #v(-20pt) 
          #image("results/csf_flow/cardiac_sas_flow/R.png", width: 100%)
        ],
        [
          #panel-label("K")
          #v(-15pt) 
          #panel-title("cardiac dispersion enhancement")
          #v(-5pt) 
          #image("results/csf_flow/cardiac_sas_flow/R_histo.png", width: 100%)
        ]
      ),

      grid(
        columns: 1, row-gutter: 5pt,
        [
          #panel-label("H")
          #v(-5pt) 
          #panel-title("respiratory-driven CSF flow")
          #v(-20pt) 
          #image("results/csf_flow/respiratory_sas_flow/csf_v.png", width: 100%)
        ],
        [
          #panel-label("L")
          #v(-15pt) 
          #panel-title("respiratory flow velocities")
          #v(-5pt)   
          #image("results/csf_flow/respiratory_sas_flow/resp_velocity_histo.png", width: 100%)
        ]
      ),

      grid(
        columns: 1, row-gutter: 5pt,
        [
          #panel-label("I")
          #v(-5pt) 
          #panel-title("respiratory dispersion enhancement")
          #v(-20pt) 
          #image("results/csf_flow/respiratory_sas_flow/R.png", width: 100%)
        ],
        [
          #panel-label("M")
          #v(-15pt) 
          #panel-title("respiratory dispersion enhancement")
          #v(-5pt) 
          #image("results/csf_flow/respiratory_sas_flow/R_histo.png", width: 100%)
        ]
      )
    ),

    // ------------------------------------------
    // BOTTOM BLOCK: Clipped Images
    // ------------------------------------------
    pad(
    right: 4cm, 
    grid(
        columns: (1fr, 1fr, 1fr, 1fr, 1fr),
        gutter: 5pt,
        [
          #panel-title("baseline")
          #v(-5pt) 
          #cropped-image("plots/modelA/modelA_overview_1-6-12-24.png", 10cm, 500%, -5.8cm, -0.3cm)
        ],
        [
          #panel-title("no advection (12 h)")
          #v(-5pt) 
          #cropped-image("plots/modelA-OnlyDispersion/modelA-OnlyDispersion_overview_1-6-12-24.png", 10cm, 500%, -5.8cm, -0.3cm)
        ],
        [
          #panel-title("low pulsatility (12 h)")
          #v(-5pt) 
          #cropped-image("plots/modelA-LowD/modelA-LowD_overview_1-6-12-24.png", 10cm, 500%, -5.8cm, -0.3cm)
        ],
        [
          #panel-title("high pulsatility (12 h)")
          #v(-5pt) 
          #cropped-image("plots/modelA-HighD/modelA-HighD_overview_1-6-12-24.png", 10cm, 500%, -5.8cm, -0.3cm)
        ],
        [
          #v(-5pt) 
          #cropped-image("plots/modelA-OnlyDispersion/modelA-OnlyDispersion_overview_1-6-12-24.png", 10cm, 500%, -11.7cm, 0.2cm)
        ]
      )
    )
  )
)