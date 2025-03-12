#set page(width: 210mm, height: 297mm, margin: 15mm)

#image("logo.png", width: 180mm)
#align(right + top)[#text(20pt, strong[RESULTS REPORT])]
#align(right + top)[*Report ID:* 2025022]

#line(length: 100%)

#block[
  #text(14pt, rgb(7, 7, 87), strong[Sample information])
]

#set text(size: 10pt)
#set table(
  stroke: none,
  align: (x, y) => (
    if x > 0 { center }
    else { left }
  )
)

#table(
  columns: 7,
  table.header(
    text(rgb(5, 5, 126))[*No.*],
    text(rgb(5, 5, 126))[*Sample ID*],
    text(rgb(5, 5, 126))[*Sample SI*],
    text(rgb(5, 5, 126))[*Gender*],
    text(rgb(5, 5, 126))[*Patient position in analysis*],
    text(rgb(5, 5, 126))[*Affection status (cause of testing)*],
    text(rgb(5, 5, 126))[*Family ID*]
  ),
  [1],
  [9-2025],
  [JD],
  [Male],
  [Probant],
  [Affected],
  [DM-152],
)

#" "
#line(length: 100%)
#block[
  #text(14pt, rgb(7, 7, 87), strong[Target information])
]

#set text(size: 8pt)
#table(
  columns: 3,
  [#table(
    align: left,
    columns: 2,
    stroke: none,
    text(rgb(5, 5, 126))[*Disease*], [*Myotonic dystrophy type 1*],
    text(rgb(5, 5, 126))[Disease abbreviation], [DM1],
    text(rgb(5, 5, 126))[OMIM ID], [#160900],
    text(rgb(5, 5, 126))[Motif complexity], [Simple (interruptions may occur)],
    text(rgb(5, 5, 126))[*Clinically relevant unit (HGVS)*], [*GCA*],
    text(rgb(5, 5, 126))[Clinically relevant unit (historical)], [CTG],
    text(rgb(5, 5, 126))[Whole motif (HGVS)], [GCA],
    text(rgb(5, 5, 126))[Whole motif (historical)], [CTG],
    text(rgb(5, 5, 126))[HGVS nomenclature (GRCh38)], [NC_000019.10:g.45770207_45770266GCA[20]],
    text(rgb(5, 5, 126))[Molecular mechanism], [Complex; Spliceopathy; ...],
    text(rgb(5, 5, 126))[Motif - Notes], [Nothing to note],
    text(rgb(5, 5, 126))[Citation (references)], [Abrakadabra et al. 2002],
  )], [#table(
        align: left,
        columns: 2,
        stroke: none,
        text(rgb(5, 5, 126))[Gene], [Dystrophia myotonica protein kinase],
        text(rgb(5, 5, 126))[Gene abbreviation], [DMPK],
        text(rgb(5, 5, 126))[Inheritance], [AD (anticipation)],
        text(rgb(5, 5, 126))[Physiological range], [5-35], 
        text(rgb(5, 5, 126))[Premutation range], [35-50],
        text(rgb(5, 5, 126))[Pathogenic range], [50+],
        text(rgb(5, 5, 126))[Grey-zone range], [-],
      )], [#table(
            align: left,
            columns: 2,
            stroke: none,
            text(rgb(5, 5, 126))[Chromosome], [19],
            text(rgb(5, 5, 126))[Gene context], [3´UTR],
            text(rgb(5, 5, 126))[Protein context], [Untranslated],
          )],
)


// #place(
//   right,
//   dx: -20pt,
//   dy: -80pt,
//   image("logo.png", width: 60mm),
// )

#show table.cell.where(y: 0): set text(
  fill: rgb(5, 5, 126),
  weight: "bold",
)

#set table(
  stroke: (x, y) => if y == 0 {
  },
  align: (x, y) => (
    if x > 0 { center }
    else { left }
  )
)

#block[
  #text(10pt, rgb(5, 5, 126), [*GRCh38 reference allele- Visualisation*])
]

#" "

#line(length: 100%)
#block[
  #text(14pt, rgb(7, 7, 87), strong[Results of analysis])
]

#table(
  columns: 7,
  table.header(
    [*Sample ID*],
    [*Sample SI*],
    [*Allele*],
    [*Repeat number prediction*],
    [*HGVS prediction*],
    [*Pathogenicity prediction*],
    [*Confidence*]
  ),
  [9-2025], [JD], [Allele 1], [5], [chr19:g.45770207_45770266GCA[5]], text(rgb(12, 152, 8))[Benign], [100.00%],
  [], [], [Allele 2], [12], [chr19:g.45770207_45770266GCA[12]], text(rgb(12, 152, 8))[Benign], [100.00%],
)

#pagebreak()

#line(length: 100%)
#block[
  #text(14pt, rgb(7, 7, 87), strong[Details of results])
]
#" "

#import "@preview/tablex:0.0.9": tablex, colspanx, rowspanx

#tablex(
  columns: (1fr, 2fr, 1fr, 1fr, 1fr, 1fr, 2fr, 1fr, 1fr, 1fr, 1fr, 2fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
  map-hlines: h => (..h, stroke: blue),  // Retaining horizontal lines in blue
  auto-vlines: false,
  fill: (_, y) => if calc.even(y) { blue.lighten(90%) },
  colspanx(2)[Sample info], colspanx(5)[Allele 1], colspanx(5)[Allele 2], colspanx(7)[Overall],

  // Rotating text headers
  rotate(-90deg, reflow: true, "Sample No."), rotate(-90deg, reflow: true,"Sample ID"),
  rotate(-90deg, reflow: true,"Repeat number"), rotate(-90deg, reflow: true,"Confidence"), rotate(-90deg, reflow: true,"Reads Spanning"), rotate(-90deg, reflow: true,"Indel errors"), rotate(-90deg, reflow: true,"Missmatch errors"),
  rotate(-90deg, reflow: true,"Repeat number"), rotate(-90deg, reflow: true,"Confidence"), rotate(-90deg, reflow: true,"Reads Spanning"), rotate(-90deg, reflow: true,"Indel errors"), rotate(-90deg, reflow: true,"Missmatch errors"),
  rotate(-90deg, reflow: true,"Confidence"), rotate(-90deg, reflow: true,"BAM - Dante - Reads"), rotate(-90deg, reflow: true,"Reads total"), rotate(-90deg, reflow: true,"Reads Spanning"), rotate(-90deg, reflow: true,"Reads Partial"),
  rotate(-90deg, reflow: true,"Indel errors"), rotate(-90deg, reflow: true,"Mismatch errors"),

  // Data rows
  [1.], [07-2025], [5], [100%], [11], [0], [0.18], [12], [100%], [5],[0],[0],[100%],[42],[21],[16],[5],[0],[0.15],
  [2.], [08-2025], [5], [100%], [11], [0], [0.18], [12], [100%], [5],[0],[0],[100%],[42],[21],[16],[5],[0],[0.15],
  [3.], [09-2025], [5], [100%], [11], [0], [0.18], [12], [100%], [5],[0],[0],[100%],[42],[21],[16],[5],[0],[0.15],
)
