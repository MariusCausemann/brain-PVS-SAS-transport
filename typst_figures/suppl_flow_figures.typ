#import "utils.typ": setup-figure, panel-label, panel-title, bg-grid
#show: setup-figure

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 20pt,
  // Column A: Cardiac CSF Velocity
  stack(
    panel-label("A"),
    panel-title("Cardiac peak flow field"),
    v(-10pt),
    image("../results/csf_flow/cardiac_sas_flow/csf_v.png", width: 113%)
  ),
  
  // Column B: Respiratory CSF Velocity
  stack(
    panel-label("B"),
    panel-title("Respiratory peak flow field"),
    v(-10pt),
    image("../results/csf_flow/respiratory_sas_flow/csf_v.png", width: 113%)
  ),
  
  // Column C: Stacked Histograms
  stack(
    spacing: 1pt,
    panel-label("C"),
    panel-title("Cardiac flow velocities"),
    image("../results/csf_flow/respiratory_sas_flow_0.5/resp_velocity_histo.png", width: 100%),
    panel-label("D"),
    panel-title("Respiratory flow velocities"),
    image("../results/csf_flow/cardiac_sas_flow_0.5/cardiac_velocity_histo.png", width: 100%)
  )
)
