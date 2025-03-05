use typst_as_lib::typst_kit_options::TypstKitFontOptions;
use typst_as_lib::TypstEngine;

use crate::analysis_family::Data as FamilyData;

pub(super) fn simple_report(data: &FamilyData) {
    let output_pdf = "typst_report.pdf";
    let typst_template = work_on_me();

    // Here be dragons
    let template = TypstEngine::builder()
        .main_file(typst_template)
        .search_fonts_with(TypstKitFontOptions::default())
        .build();

    let doc = template.compile().output
        .expect("typst::compile() returned an error!");

    let options = Default::default();

    let pdf = typst_pdf::pdf(&doc, &options).expect("Could not generate pdf.");
    std::fs::write(output_pdf, pdf).expect("Could not write pdf.");
}

fn work_on_me() -> String {
    let typst_template: &str = "
== Hello World

In this report, we will explore the
various factors that influence fluid
dynamics in glaciers and how they
contribute to the formation and
behaviour of these natural structures.

+ The climate
  - Temperature
  - Precipitation
+ The topography
+ The geology
    ";
    return typst_template.to_string();
}
