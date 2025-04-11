#set page(
  paper: "a4",
  margin: (top: 25mm, right: 10mm, bottom: 25mm, left: 10mm)
)

#let dark_blue = rgb(7, 7, 87)
#let blue = rgb(5, 5, 126)
#let light_blue = rgb(195, 215, 255)
#let r = 2mm

#let benign = rect(radius: 50%, fill: green)[Benign]
#let pathogenic = rect(radius: 50%, fill: red)[Pathogenic]
#let premutation = rect(radius: 50%, fill: yellow)[Premutation]

// some useful functions for debugging and prototyping
#let dbox(content) = box(stroke: red)[content] /* this is useful for debuging layouts */
#let placeholder(width, height) = box(
  width: width, height: height, fill: gray, stroke: gray
)[#align(center + horizon)[Image]]

// ----------------------------------------------------------------------------

#table(
  columns: (60%, 40%),
  stroke: none,
  inset: 0mm,
  [#placeholder(100%, 10%)],
  [
    #align(right)[#text(20pt)[*ONE-PAGE REPORT*]]
    #align(right)[*Report ID:* {{ report_id }}]
  ]
)

#v(5mm)
#rect(width: 100%,fill: light_blue, radius: r, inset: 2.5mm)[

  #rect(radius: r, fill: white)[
    #table(columns: (50%, 50%), inset: 1mm, fill: white, stroke: none)[
      #align(left)[
        #text(16pt, dark_blue)[*Proband*]
        #h(0.5cm)
        #text(15pt)[{{ proband_id }}]
      ]
    ][
      #align(right)[
        #text(16pt, dark_blue)[*Family ID*]
        #h(0.5cm)
        #text(14pt)[{{ family_id }}]
      ]
    ]
  ]

  #rect(fill: white, radius: r, height: 1cm)[
    #text(14pt, dark_blue)[*Family structure*]
  ]

  #v(-8mm) /* 5mm seems to be the diff between sections */
  #rect(fill: white, radius: r)[
    #table(
      columns: (3fr, 9fr, 9fr, 6fr, 8fr, 16fr, 13fr), /* number of letters in header seems like a good heuristic */
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
        text(blue)[*Date of birth*]
      ),
      [{{ row[0] }}],
      [{{ row[1] }}],
      [{{ row[2] }}],
      [{{ row[3] }}],
      [{{ row[4] }}],
      [{{ row[5] }}],
      [{{ row[6] }}],
    )
  ]

  #rect(fill: white, radius: r, height: 1cm)[
    #text(14pt, dark_blue)[*Case description*]
  ]

  #v(-8mm) /* 5mm seems to be the diff between sections */
  #rect(fill: white, radius: r)[
    // This is some serious gourmet feature
    #show table.cell.where(x: 0): set text(
      fill: blue,
      weight: "bold",
    )

    #table(
      columns: 2,
      stroke: none,
      align: left,
      [Suspected diagnosis],              [{{ sus_diag }}],
      [Requesting person],                [{{ req_person }}],
      [Reason for testing],               [{{ reason }}],
      [Requesting facility],              [{{ req_facility }}],
      [Other health conditions],          [{{ health_cond }}],
      [Phenotype in HPO terms],           [{{ hpo_terms }}],
      [Requested target(s) of analysis],  [{{ req_targets }}],
      [Patient note],                     [{{ note }}],
    )
  ]
  
  #rect(fill: white, radius: r, height: 1cm)[
    #text(14pt, dark_blue)[*Analysis resume*]
  ]
  #v(-8mm)
  #rect(fill: white, radius: r, width: 100%)[
    #table(
      columns: 6,
      stroke: none,
      [#text(blue)[*No. of requested targets:*]], [{{ n_req_targets }}],
      [#text(blue)[*No. of QC passed loci:*]], [{{ n_loci_qc_pass }}],
      [#text(blue)[*No. of QC failed loci:*]], [{{ n_loci_qc_fail }}],
    )
    #v(-4.7mm)
    #table(
      columns: 2,
      stroke: none,
      [#text(blue)[*Failure reasons:*]], [{{ fail_reason }}]
    )
  ]

  #rect(fill: white, radius: r, height: 1cm)[
    #text(14pt, dark_blue)[*Results interpretation*]
  ]
  #v(-8mm)
  #rect(fill: white, radius: r, width: 100%)[
    {{ interpretation }}
    #v(2mm)
  ]

  #rect(fill: white, radius: r, height: 1cm)[
    #text(14pt, dark_blue)[*Clinically relevant findings*]
  ]
  #v(-8mm)
  #rect(fill: white, radius: r, width: 100%)[
    #table(
      columns: 7,
      stroke: none,
      align: left,
      table.header(
        [#text(blue)[*No.*]],
        [#text(blue)[*Sample ID*]],
        [#text(blue)[*Sample SI*]],
        [#text(blue)[*Allele*]],
        [#text(blue)[*No. of repeats*]], /* [Repeat number \ (revision)], */
        [#text(blue)[*Pathogenicity*]], /* [HGVS nomenclature \ (revision)], */
        [#text(blue)[*HGVS nomenclature*]], /* [Pathogenicity \ (revision)], */
      ),
        [{{ row[0] }}], [{{ row[1] }}], [{{ row[2] }}],
        [Allele 1], [{{ a1_repnum }}], [{{ a1_status }}], [{{ a1_nom }}],
        [], [], [],
        [Allele 2], [{{ a2_repnum }}], [{{ a2_status }}], [{{ a2_nom }}],
    )
  ]
]
