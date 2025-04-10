#set page(paper: "a4", margin: 15mm)

#let dark_blue = rgb(7, 7, 87)
#let blue = rgb(5, 5, 126)
#let light_blue = rgb(195, 215, 255)
#let green = rgb(178, 255, 108)
#let r = 1mm

#align(right + top)[#text(20pt, strong[ONE-PAGE REPORT])]
#align(right + top)[*Report ID:* {2025022}]

#place(right, dy: -8pt, dx: 9pt, rect(width: 190mm, height: 80%, fill: light_blue, radius: r))

#v(0.2cm)
#place(right, dy: -6pt, dx: 2pt, rect(width: 185mm, height: 25pt, fill: white, radius: r))
#block[
  #text(16pt, dark_blue, strong[Proband]) #h(1cm)
  #text(15pt)[{9-2025}] #h(7cm)
  #text(16pt, dark_blue, strong[Family ID]) #h(1cm)
  #text(14pt)[{DM-152}]
]

#v(0.4cm)
#place(right, dy: -5pt, dx: -389pt, rect(width: 47mm, height: 60pt, fill: white, radius: r))
#place(right, dy: 15pt, dx: 2pt, rect(width: 185mm, height: 60pt, fill: white, radius: r))
#block[#text(14pt, dark_blue, strong[Family structure])]

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
    text(blue)[*No.*],
    text(blue)[*Sample ID*],
    text(blue)[*Sample SI*],
    text(blue)[*Gender*],
    text(blue)[*Patient position in analysis*],
    text(blue)[*Affection status \ (cause of testing)*],
    text(blue)[*Date of birth*]
  ),
  [{1}],
  [{9-2025}],
  [{JD}],
  [{Male}],
  [{Proband}],
  [{Affected}],
  [{1954-12-02}],
)

#v(0.4cm)
#place(right, dy: -5pt, dx: -389pt, rect(width: 47mm, height: 60pt, fill: white, radius: r))
#place(right, dy: 15pt, dx: 2pt, rect(width: 185mm, height: 175pt, fill: white, radius: r))
#block[
  #text(14pt, dark_blue, strong[Case description])
]

#show table.cell.where(x: 0): set text(
  fill: rgb(5, 5, 126),
  weight: "bold",
)

#show table.cell.where(x: 2): set text(
  fill: rgb(5, 5, 126),
  weight: "bold",
)

#set table(
  stroke: (x, y) => if y == 0 {},
  align: (x, y) => (
    if x > 0 { center }
    else { left }
  )
)

#table(
  columns: 4,
  align: left,
  [Suspected diagnosis], [Myotonic dystrophy], [Requesting person], [Dr. Umberto Eco],
  
  [Reason for testing],
  [Patient has long term muscle pain (proximal muscles), cataract surgery 5 years ago, etc.],
  [Requesting facility],
  [Milano Faculty Hospital],
  
  [Other health conditions], [High blood pressure], [Phenotype in HPO terms], [proximal muscle pain, cataract],
  [Requested target(s) of analysis], [DMPK (DM1); CNBP (DM2)], [], [],
  [Patient note], [Nothing to note], [], []
)

#v(0.8cm)
#place(right, dy: -5pt, dx: -389pt, rect(width: 47mm, height: 60pt, fill: white, radius: r))
#place(right, dy: 15pt, dx: 2pt, rect(width: 185mm, height: 45pt, fill: white, radius: r))
#block[
  #text(14pt, rgb(7, 7, 87), strong[Analysis resume])
]

#show table.cell.where(x: 4): set text(
  fill: rgb(5, 5, 126),
  weight: "bold",
)

#table(
  columns: 6,
  [No. of requested target(s):], [2], [No. of QC passed loci], [2], [No. of QC failed loci], [0],
  [Failure reason(s)], [None], [], [], [], []
)

#v(0.4cm)
#place(right, dy: -5pt, dx: -355pt, rect(width: 59mm, height: 30pt, fill: white, radius: r))
#place(right, dy: 15pt, dx: 2pt, rect(width: 185mm, height: 35pt, fill: white, radius: r))
#block[
  #text(14pt, rgb(7, 7, 87), strong[Results interpretation])
]
#block[From the 2 analyzed target loci all 2 passed the QC filter and all 2 were interpretable. We identified no pathogenic repeat structures in these loci... However, the repeat structure in the DM1 (DMPK) CTG motif was found to be atypical...]

#v(0.5cm)
#place(right, dy: -5pt, dx: -324pt, rect(width: 70mm, height: 60pt, fill: white, radius: r))
#place(right, dy: 15pt, dx: 2pt, rect(width: 185mm, height: 90pt, fill: white, radius: r))
#block[
  #text(14pt, rgb(7, 7, 87), strong[Clinically relevant findings])
]

#show table.cell.where(x: 2): set text(
  fill: black,
  weight: "regular",
)

#show table.cell.where(x: 4): set text(
  fill: black,
  weight: "regular",
)

#show table.cell.where(y: 0): set text(
  fill: rgb(5, 5, 126),
  weight: "bold",
)

#table(
  columns: 7,
  table.header(
    [No.], [Sample ID], [Sample SI], [Allele], [Repeat number (revision)], [HGVS nomenclature (revision)], [Pathogenicity (revision)],
  ),
    [{1}], [{9-2025}], [{JD}], [Allele 1], [5],  [chr19:g.45770207_45770266GCA[5]],  rect(radius: 50%, fill: green)[Benign],
    [],  [],       [],   [Allele 2], [12], [chr19:g.45770207_45770266GCA[12]], rect(radius: 50%, fill: green)[Benign],
)
