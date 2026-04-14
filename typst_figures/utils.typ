// theme.typ

#let setup-figure(body) = {
  // Global settings for the figure
  set page(width: 21cm, height: auto, margin: 1cm)
  set text(font: "Nimbus Sans", size: 8pt)
  
  body
}

#let panel-label(content, nudge-dx: 0pt, nudge-dy: 0pt) = align(top + left)[
  #move(dx: nudge-dx, dy: nudge-dy)[#text(weight: "bold", size: 12pt)[#content]]
]

#let panel-title(content) = align(center)[
  #text(weight: "regular", size: 10pt)[#content]
]

#let inline-header(lbl, ttl) = align(left)[
  #box[#text(weight: "bold", size: 12pt)[#lbl]]
  #h(4pt)
  #box[#text(weight: "regular", size: 10pt)[#ttl]]
  #v(2pt) // slight buffer before the image
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

#let y-label(content, nudge: 0cm) = align(center + horizon)[
  #move(dy: nudge)[
    #rotate(-90deg, reflow: true)[#text(size: 6.5pt)[#content]]
  ]
]