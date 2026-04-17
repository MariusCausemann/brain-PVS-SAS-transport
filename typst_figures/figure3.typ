#import "utils.typ": setup-figure, panel-label, panel-title, inline-header, cropped-image
#show: setup-figure


#let move_x = -7.2cm

// ==========================================
// MAIN FIGURE LAYOUT
// ==========================================
#figure(
  grid(
    columns: 1,
    row-gutter: 5pt, 
    // ------------------------------------------
    // TOP BLOCK: A/B | C/D | E/F
    // ------------------------------------------
    grid(
      columns: (1fr, 1fr, 1fr),
      gutter: 0pt,
      
      // Column 1: Production-driven
      grid(
        columns: 1, row-gutter: 5pt,
        [
          #panel-label("A")
          #v(-15pt) 
          #panel-title("Production-driven CSF flow")
          #v(-25pt) 
          #image("../results/csf_flow/sas_flow/csf_v.png", width: 120%)
        ],
        [
          #panel-label("B")
          #v(-15pt) 
          #panel-title("Production-driven CSF velocities")
          #v(-5pt)
          #image("../results/csf_flow/sas_flow/prod_velocity_histo.png", width: 100%)
        ]
      ),

      // Column 2: Cardiac Dispersion
      grid(
        columns: 1, row-gutter: 5pt,
        [
          #panel-label("C")
          #v(-15pt) 
          #panel-title("Cardiac dispersion enhancement")
          #v(-25pt) 
          #image("../results/csf_flow/cardiac_sas_flow/R.png", width: 120%)
        ],
        [
          #panel-label("D")
          #v(-15pt) 
          #panel-title("Cardiac dispersion enhancement")
          #v(-5pt) 
          #image("../results/csf_flow/cardiac_sas_flow/R_histo.png", width: 100%)
        ]
      ),

      // Column 3: Respiratory Dispersion
      grid(
        columns: 1, row-gutter: 5pt,
        [
          #panel-label("E")
          #v(-15pt) 
          #panel-title("Respiratory dispersion enhancement")
          #v(-25pt) 
          #image("../results/csf_flow/respiratory_sas_flow/R.png", width: 120%)
        ],
        [
          #panel-label("F")
          #v(-15pt) 
          #panel-title("Respiratory dispersion enhancement")
          #v(-5pt) 
          #image("../results/csf_flow/respiratory_sas_flow/R_histo.png", width: 100%)
        ]
      )
    ),

    // ------------------------------------------
    // BOTTOM BLOCK: Clipped Images (G, H, I, J)
    // ------------------------------------------
    pad(
      right: 0.5cm, 
      left : 0.0cm,
      grid(
        columns: (1fr, 1fr, 1fr, 1fr, 1fr),
        gutter: 5pt,
        [
          #inline-header("G", "Baseline (12 h)")
          #v(-10pt)
          #cropped-image("../plots/modelA/modelA_overview_1-6-12-24.png", 11.2cm, 500%, move_x, -0.3cm)
        ],
        [
          #inline-header("H", "No advection")
          #v(-10pt)
          #cropped-image("../plots/modelA-OnlyDispersion/modelA-OnlyDispersion_overview_1-6-12-24.png", 11.2cm, 500%, move_x, -0.3cm)
        ],
        [
          #inline-header("I", "Low pulsatility")
          #v(-10pt)
          #cropped-image("../plots/modelA-LowD/modelA-LowD_overview_1-6-12-24.png", 11.2cm, 500%, move_x, -0.3cm)
        ],
        [
          #inline-header("J", "High pulsatility")
          #v(-10pt)
          #cropped-image("../plots/modelA-HighD/modelA-HighD_overview_1-6-12-24.png", 11.2cm, 500%, move_x, -0.3cm)
        ],
        [
          #inline-header("", "")
          #v(-10pt)
          #cropped-image("../plots/modelA-OnlyDispersion/modelA-OnlyDispersion_overview_1-6-12-24.png", 11.2cm, 500%, -14.5cm, -0.3cm)
        ]
      )
    ),
    v(-25pt)
  )
)