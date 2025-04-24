#set page(margin: 10mm)
// This is correct font to use
#set text(font: "Fira Sans", stretch: 50%, size: 8pt)

#let dark_blue = rgb(7, 7, 87)
#let blue = rgb(5, 5, 126)
#let light_blue = rgb(195, 215, 255)
#let r = 1mm

#let benign = rect(radius: 50%, fill: rgb(166, 227, 161))[Benign]
#let premutation = rect(radius: 50%, fill: rgb(249, 226, 175))[Premutation]
#let pathogenic = rect(radius: 50%, fill: rgb(243, 139, 168))[Pathogenic]
#let unknown = rect(radius: 50%, fill: rgb(205, 214, 244))[Unknown]

// Oklab interpolation between benign/pathogenic and unknown
#let likely_benign = rect(radius: 50%, fill: rgb(185, 221, 204))[Likely Benign]
#let likely_pathogenic = rect(radius: 50%, fill: rgb(224, 178, 205))[Likely Pathogenic]

#let positive(x) = rect(radius: 50%, fill: rgb(166, 227, 161), x)
#let negative(x) = rect(radius: 50%, fill: rgb(243, 139, 168), x)

// some useful functions for debugging and prototyping
#let dbox(content) = box(stroke: red)[content] /* this is useful for debuging layouts */

// Determine if a path is a valid image (using magic) and either show it or use some placeholder.
#let maybe_image(path, width, height) = context {
  let context_function = (context { }).func()
  let first_time = query(context_function).len() == 0
  let path_label = label(path)
  let used_path = query(path_label).len() > 0
  if first_time or used_path {
    [#image(path, height: height)#path_label]
  } else {
    [#box(width: width, height: height, fill: gray)[#align(center + horizon, path)]]
  }
}

// ---------------------------------------------------------------------------------

#table(
  columns: (60%, 40%),
  stroke: none,
  inset: 0mm,
  [#maybe_image("logo", 100%, 10%)],
  [
    #align(right)[#text(20pt)[*RESULTS REPORT*]]
    #align(right)[*Report ID:* {{g.report_id}}]
  ]
)

#rect(width: 100%,fill: light_blue, radius: r, inset: 2.5mm)[
  #rect(fill: white, radius: r, height: 1cm)[#text(14pt, dark_blue)[*Sample information*]]
  #v(-7mm)
  #rect(fill: white, radius: r, width: 100%)[
    #table(
      columns: (3fr, 9fr, 9fr, 6fr, 8fr, 16fr, 13fr),
      stroke: none,
      align: (x, _) => (
        if x > 0 { center }
        else { left }
      ),
      table.header(
        text(blue)[*No.*],
        text(blue)[*Sample ID*],
        text(blue)[*Sample SI*],
        text(blue)[*Gender*],
        text(blue)[*Position*], /* [*Patient position in analysis*] */
        text(blue)[*Affection status*], /* [*Affection status \ (cause of testing)*] */
        text(blue)[*Family ID*]
      ),
      [1.],
      [{{m.pid}}],
      [{{m.sid}}],
      [{{m.gender}}],
      [Proband],
      [{{m.status}}],
      [{{m.fid}}],
    )
  ]
]

#rect(width: 100%,fill: light_blue, radius: r, inset: 2.5mm)[
  #rect(fill: white, radius: r, height: 1cm)[#text(14pt, dark_blue)[*Target information*]]
  #v(-7mm)
  #rect(fill: white, radius: r, width: 100%)[
    #[
      #show table.cell: it => {
        if it.x == 0 or it.x == 2 or it.x == 4 { set text(blue); strong(it) }
        else { it }
      }
      #table(
        columns: 6,
        stroke: none,
        [Disease], [{{h.name}}], [Gene], [{{h.gene}}], [Chromosome], [{{h.chr}}],
        [Disease abbreviation], [{{h.abbr}}], [Gene abbreviation], [{{h.gene_abbr}}], [Gene context], [{{h.gene_ctx}}],
        [OMIM ID], [\{{h.omim_id}}], [Inheritance], [{{h.inheritance}}], [Protein context], [{{h.prot_ctx}}],
        [Motif complexity], [{{h.motif_cpx}}],
      )
    ]
    #table(
      columns: 2,
      stroke: none,
      align: left + horizon,
      [
        #show table.cell: it => {
          if it.x == 0 or it.x == 2 { set text(blue); strong(it) }
          else { it }
        }
        #table(
          columns: 4,
          stroke: none,
          [Module], [{{h.module}}],
          [Clinically relevant unit (HGVS)], [{{h.unit_hgvs}}],
          [Physiological range], [{{h.physiological}}],
          [Clinically relevant unit (historical)], [{{h.unit_hist}}],
          [Premutation range], [{{h.premutation}}],
          [Whole motif (HGVS)], [{{h.motif_hgvs}}],
          [Pathogenic range], [{{h.pathogenic}}],
          [Whole motif (historical)], [{{h.motif_hist}}],
        )
      ],
      [#maybe_image("{{h.dist_image}}", 100%, 4cm)]
    )
    #[
      #show table.cell: it => {
          if it.x == 0 { set text(blue); strong(it) }
          else { it }
      }
      #table(
        columns: 2,
        stroke: none,
        align: horizon,
        [GRCh38 reference allele \ (HGVS nomenclature)], [{{h.ref_allele_hgvs}}],
        [GRCh38 reference allele \ (Visualization)], [#maybe_image("visualization", 100%, 12pt)],
        [Molecular mechanism], [{{h.mechanism}}],
        [Motif - Notes], [{{h.notes}}],
        [Citation (references)], [{{h.citations}}]
      )
    ]
  ]
]

#rect(width: 100%,fill: light_blue, radius: r, inset: 2.5mm)[
  #rect(fill: white, radius: r, height: 1cm)[#text(14pt, dark_blue)[*Results of analysis*]]
  #v(-7mm)
  #rect(fill: white, radius: r, width: 100%)[
    #table(
      columns: (7mm, auto, auto, auto, auto, auto, 1fr, auto), /* 8, */ 
      stroke: none,
      align: horizon,
      inset: (x: 2mm, y: 0.75mm),
      table.header(
        text(blue)[*No.*],
        text(blue)[*Module*],
        text(blue)[*Sample ID*],
        text(blue)[*Sample SI*],
        text(blue)[*Allele*],
        text(blue)[*No. repeats \ (revision)*]/*[*Repeat number \ (revision)*]*/,
        text(blue)[*HGVS nomenclature \ (revision)*],
        text(blue)[*Pathogenicity \ (revision)*],
      ),
      [{G-1.}], [{G-DM2-0}], [{{m.pid}}], [{{m.sid}}], [Allele 1], [{R-5}],  [{R-chr19:g.45770207_45770266GCA[5]}], benign,
      [],       [],          [],           [],       [Allele 2], [{R-12}], [{R-chr19:g.45770207_45770266GCA[12]}], likely_benign,
      [{G-2.}], [{G-DM2-1}], [{{m.pid}}], [{{m.sid}}], [Allele 1], [{R-5}],  [{R-chr19:g.45770207_45770266GCA[5]}], premutation,
      [],       [],          [],           [],       [Allele 2], [{R-12}], [{R-chr19:g.45770207_45770266GCA[12]}], likely_pathogenic,
      [{G-3.}], [{G-DM2-2}], [{{m.pid}}], [{{m.sid}}], [Allele 1], [{R-5}],  [{R-chr19:g.45770207_45770266GCA[5]}], pathogenic,
      [],       [],          [],           [],       [Allele 2], [{R-12}], [{R-chr19:g.45770207_45770266GCA[12]}], unknown,
    )
  ]
  #rect(fill: white, radius: r, width: 100%)[
    #table(
      columns: (20mm, 1fr, auto, 1fr),
      stroke: none,
      align: horizon,
      inset: (x: 2mm, y: 0.5mm),
      [#text(blue)[*Results QC*]], [{R-#positive[Passed]}],
      [#text(blue)[*Included in One-page report?*]], [{R-#negative[No]}],
    )
    #v(-2mm)
    #box(width: 20mm, inset: (x: 2mm, y: 0mm))[#text(blue)[*Failure reason*]]
    #box(inset: (x: 2mm, y: 0mm))[
      {R-Low read number and some veeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeery long reason.}
    ]
  ]
]

// #pagebreak()

#let h1(x) = rotate(-90deg, reflow: true)[#text(blue, x)];
#let mc(x) = table.cell(rowspan: 2)[#x];

#rect(width: 100%,fill: light_blue, radius: r, inset: 2.5mm)[
  #rect(fill: white, radius: r, height: 1cm)[#text(14pt, dark_blue)[*Details of results*]]
  #v(-7mm)
  #rect(fill: white, radius: r, width: 100%)[
    #align(center)[#text(10pt, blue)[*Dante predictions*]]
    #table(
      columns: (
        auto, auto, auto, auto, 1fr,
        auto, auto, auto, auto, auto, auto, 1fr,
        auto, auto, auto, auto, auto, auto, auto, 1fr,
        auto, auto
      ),
      align: (x, y) => {
        if y == 1 { center + bottom }
        else { center + horizon }
      },
      stroke: none,
      inset: (x: 2mm, y: 0.75mm),
      table.cell(colspan: 4)[#text(blue)[*Sample info*]],
      [],
      table.cell(colspan: 6)[#text(blue)[*Results*]],
      [],
      table.cell(colspan: 7)[#text(blue)[*Overall statistics*]],
      [],
      table.cell(colspan: 2)[#text(blue)[*Flags*]],

      h1[*No.*],
      h1[*Module*],
      h1[*Sample ID*],
      h1[*Sample SI*],
      [],
      h1[*Allele*],
      h1[*No. repeats*],
      h1[*Pathogenicity*],
      h1[*Reads spanning*],
      h1[*Indel errors*],
      h1[*Mismatch errors*],
      [],
      h1[*Confidence*],
      h1[*BAM-Dante - Reads*],
      h1[*Reads total*],
      h1[*Reads spanning*],
      h1[*Reads partial*],
      h1[*Indel errors*],
      h1[*Mismatch errors*],
      [],
      h1[*Alert of edited results*],
      h1[*Results QC*],

      table.hline(y: 2, stroke: 0.1mm),
      mc[{G-1.}], mc[{D-DM2-0}], mc[{{m.pid}}], mc[{{m.sid}}], mc[],
      [Allele 1], [{D-5}],  [{?-#benign}],        [{D-10}], [{D-2}], [{D-1}],
      mc[], mc[{D-90%}], mc[{D-10}], mc[{D-10}], mc[{D-7}], mc[{D-3}], mc[{D-0.3}], mc[{D-0.5}], mc[], mc[{R-!}], mc[{R-OK}],
      [Allele 2], [{D-12}], [{?-#likely_benign}], [{D-2}],  [{D-2}], [{D-1}],

      table.hline(y: 4, stroke: 0.1mm),
      mc[2.], mc[DM2-1], mc[{{m.pid}}], mc[{{m.sid}}], mc[],
      [Allele 1], [8],  [#pathogenic],        [10], [2], [1],
      mc[], mc[90%], mc[10], mc[10], mc[7], mc[3], mc[0.3], mc[0.5], mc[], mc[!], mc[OK],
      [Allele 2], [20], [#unknown], [2],  [2], [1],

      table.hline(y: 6, stroke: 0.1mm),
      mc[2.], mc[DM2-2], mc[{{m.pid}}], mc[{{m.sid}}], mc[],
      [Allele 1], [17],  [#benign],        [10], [2], [1],
      mc[], mc[90%], mc[10], mc[10], mc[7], mc[3], mc[0.3], mc[0.5], mc[], mc[!], mc[OK],
      [Allele 2], [17], [#benign], [2],  [2], [1],

      table.hline(y: 8, stroke: 0.1mm),
    )

    #v(1em)
    #align(center)[#text(10pt, blue)[*Revised predictions*]]
    #table(
      columns: (
        auto, auto, auto, auto, 1fr,
        auto, auto, auto, auto, auto, auto, 1fr,
        auto, auto, auto, auto, auto, auto, auto, 1fr,
        auto, auto
      ),
      align: (x, y) => {
        if y == 1 { center + bottom }
        else { center + horizon }
      },
      stroke: none,
      inset: (x: 2mm, y: 0.75mm),
      table.cell(colspan: 4)[#text(blue)[*Sample info*]],
      [],
      table.cell(colspan: 6)[#text(blue)[*Results*]],
      [],
      table.cell(colspan: 7)[#text(blue)[*Overall statistics*]],
      [],
      table.cell(colspan: 2)[#text(blue)[*Flags*]],

      h1[*No.*],
      h1[*Module*],
      h1[*Sample ID*],
      h1[*Sample SI*],
      [],
      h1[*Allele*],
      h1[*No. repeats*],
      h1[*Pathogenicity*],
      h1[*Reads spanning*],
      h1[*Indel errors*],
      h1[*Mismatch errors*],
      [],
      h1[*Confidence*],
      h1[*BAM-Dante - Reads*],
      h1[*Reads total*],
      h1[*Reads spanning*],
      h1[*Reads partial*],
      h1[*Indel errors*],
      h1[*Mismatch errors*],
      [],
      h1[*Alert of edited results*],
      h1[*Results QC*],

      table.hline(y: 2, stroke: 0.1mm),
      mc[1.], mc[DM2-0], mc[{{m.pid}}], mc[{{m.sid}}], mc[],
      [Allele 1], [5],  [#benign],        [10], [2], [1],
      mc[], mc[90%], mc[10], mc[10], mc[7], mc[3], mc[0.3], mc[0.5], mc[], mc[!], mc[OK],
      [Allele 2], [12], [#likely_benign], [2],  [2], [1],

      table.hline(y: 4, stroke: 0.1mm),
      mc[2.], mc[DM2-1], mc[{{m.pid}}], mc[{{m.sid}}], mc[],
      [Allele 1], [8],  [#pathogenic],        [10], [2], [1],
      mc[], mc[90%], mc[10], mc[10], mc[7], mc[3], mc[0.3], mc[0.5], mc[], mc[!], mc[OK],
      [Allele 2], [20], [#unknown], [2],  [2], [1],

      table.hline(y: 6, stroke: 0.1mm),
      mc[2.], mc[DM2-2], mc[{{m.pid}}], mc[{{m.sid}}], mc[],
      [Allele 1], [17],  [#benign],        [10], [2], [1],
      mc[], mc[90%], mc[10], mc[10], mc[7], mc[3], mc[0.3], mc[0.5], mc[], mc[!], mc[OK],
      [Allele 2], [17], [#benign], [2],  [2], [1],

      table.hline(y: 8, stroke: 0.1mm), 
    )
  ]
]

#v(-5mm)
#rect(width: 100%,fill: light_blue, radius: r, inset: 2.5mm)[
  #rect(fill: white, radius: r, width: 100%)[
    #table(
      columns: (1fr, 1fr, 50%, 1fr),
      stroke: none,
      // align: center + horizon,
      align: (horizon, horizon, center + horizon, horizon),
      [], /* [#text(blue)[*Family tree*]],*/
      [#text(blue)[*Sample ID \ Genotype prediction \ Genotype revision*]],
      [#text(blue)[*Histogram*]],
      [#text(blue)[*Nomenclatures*]]
    )
  ]
]
#v(-5mm)
#rect(width: 100%,fill: light_blue, radius: r, inset: 2.5mm)[
  #rect(fill: white, radius: r, width: 100%)[
    #table(
      columns: (1fr, 1fr, 50%, 1fr),
      stroke: none,
      align: horizon,
      [],
      [{{m.pid}} \ 5/12 \ 5/12],
      [#maybe_image("DM2_0_histogram.png", 100%, 65mm)],
      [10x GCA[10] \ 7x GCA[4]]
    )
  ]
]
#v(-5mm)
#rect(width: 100%,fill: light_blue, radius: r, inset: 2.5mm)[
  #rect(fill: white, radius: r, width: 100%)[
    #table(
      columns: (1fr, 1fr, 50%, 1fr),
      stroke: none,
      align: horizon,
      [],
      [{{m.pid}} \ 5/12 \ 5/12],
      [#maybe_image("DM2_1_histogram.png", 100%, 65mm)],
      [10x GCA[10] \ 7x GCA[4]]
    )
  ]
]
#v(-5mm)
#rect(width: 100%,fill: light_blue, radius: r, inset: 2.5mm)[
  #rect(fill: white, radius: r, width: 100%)[
    #table(
      columns: (1fr, 1fr, 50%, 1fr),
      stroke: none,
      align: horizon,
      [],
      [{{m.pid}} \ 5/12 \ 5/12],
      [#maybe_image("DM2_2_histogram.png", 100%, 65mm)],
      [10x GCA[10] \ 7x GCA[4]]
    )
  ]
]
