// theme.typ

#let setup-figure(body) = {
  // Global settings for the figure
  set page(width: 21cm, height: auto, margin: 1cm)
  set text(font: "Nimbus Sans", size: 8pt)
  // Temporary alignment grid
  body
}

#let bg-grid() = {
  set page(background: {
    place(top + left)[
      #grid(
        columns: (1cm,) * 30, 
        rows: (1cm,) * 60,
        stroke: 0.5pt + luma(200) // Light gray lines
      )
    ]
  })
}

#let panel-label(content, nudge-dx: 0pt, nudge-dy: 0pt) = align(top + left)[
  #move(dx: nudge-dx, dy: nudge-dy)[#text(weight: "bold", size: 12pt)[#content]]
]

#let panel-title(content, fontsize: 9pt, alignment: center, dx: 0pt) = {
  align(alignment)[
    #move(dx: dx)[
      #text(weight: "regular", size: fontsize)[#content]
    ]
  ]
}

#let inline-header(lbl, ttl) = align(left)[
  #box[#text(weight: "bold", size: 12pt)[#lbl]]
  #h(4pt)
  #box[#text(weight: "regular", size: 9pt)[#ttl]]
  #v(2pt) // slight buffer before the image
]


#let cropped-image(filename, box-height, zoom, move-x, move-y, width: 100%) = {
  box(
    width: width, 
    height: box-height, 
    clip: true,
    place(
      dx: move-x, 
      dy: move-y, 
      image(filename, width: zoom)
    )
  )
}

#let crop-image(path, t, r, b, l, bwidth: 100%, ..args) = {
  box(
    clip: true,
    width : bwidth,
    inset: (top: -t, right: -r, bottom: -b, left: -l),
    image(path, ..args)
  )
}

#let y-label(content, nudge: 0cm) = align(center + horizon)[
  #move(dy: nudge)[
    #rotate(-90deg, reflow: true)[#text(size: 6.5pt)[#content]]
  ]
]

#let x-label(content, nudge: 0cm) = align(center + horizon)[
  #move(dy: nudge)[
    #text(size: 6.5pt)[#content]
  ]
]