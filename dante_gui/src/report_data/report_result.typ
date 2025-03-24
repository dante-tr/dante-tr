#set page(margin: 15mm) // A4 is default

#image("logo.png", width: 180mm)
#align(right + top)[#text(20pt, strong[RESULTS REPORT])]
#align(right + top)[*Report ID:* 2025022]

#place(right, dy: -8pt, dx: 9pt, rect(width: 190mm, height: 86pt, fill: rgb(195, 215, 255), radius: 15%))
#place(right, dy: 20pt, rect(width: 182mm, height: 51pt, fill: rgb(255, 255, 255), radius: 15%))
#block[
  #text(14pt, rgb(7, 7, 87), strong[Sample information])]

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
    text(rgb(5, 5, 126))[*Affection status*],
    text(rgb(5, 5, 126))[*Family ID*]
  ),
  [1],[7-2025],[JD],[Male],[Probant],[Affected],[DM-152],
  [2],[7-2025],[JD],[Male],[Probant],[Affected],[DM-152],
)

#place(right, dy: 10pt, dx: 9pt, rect(width: 190mm, height: 320pt, fill: rgb(195, 215, 255), radius: 5%))
#place(right, dy: 38pt, dx: -2pt, rect(width: 182mm, height: 285pt, fill: rgb(255, 255, 255), radius:5%))
#" "
#block[#text(14pt, rgb(7, 7, 87), strong[Target information])]
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
    text(rgb(5, 5, 126))[*Clinically relevant unit (HGVS)*], [**],
    text(rgb(5, 5, 126))[Clinically relevant unit (historical)], [],
    text(rgb(5, 5, 126))[Whole motif (HGVS)], [],
    text(rgb(5, 5, 126))[Whole motif (historical)], [],
    text(rgb(5, 5, 126))[HGVS nomenclature (GRCh38)], [],
    text(rgb(5, 5, 126))[Molecular mechanism], [Complex; Spliceopathy;],
    text(rgb(5, 5, 126))[Motif - Notes], [Nothing to note],
    text(rgb(5, 5, 126))[Citation (references)], [Abrakadabra et al. 2002],
  )], [#table(
    align: left,
    columns: 2,
    stroke: none,
    text(rgb(5, 5, 126))[Gene], [Dystrophia myotonica protein kinase],
    text(rgb(5, 5, 126))[Gene abbreviation], [DMPK],
    text(rgb(5, 5, 126))[Inheritance], [AD (anticipation)],
    text(rgb(5, 5, 126))[Physiological range], [5-30],
    text(rgb(5, 5, 126))[Premutation range], [36-50],
    text(rgb(5, 5, 126))[Pathogenic range], [51+],
    text(rgb(5, 5, 126))[Grey-zone range], [0-4,31-35],
  )], [#table(
    align: left,
    columns: 2,
    stroke: none,
    text(rgb(5, 5, 126))[Chromosome], [],
    text(rgb(5, 5, 126))[Gene context], [3´UTR],
    text(rgb(5, 5, 126))[Protein context], [Untranslated],
  )],
)
