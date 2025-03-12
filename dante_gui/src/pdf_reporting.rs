use typst_as_lib::typst_kit_options::TypstKitFontOptions;
use typst_as_lib::TypstEngine;
use typst_as_lib::package_resolver::PackageResolver;
use typst_as_lib::package_resolver::FileSystemCache;

use std::path::PathBuf;

use crate::analysis_family::Data as FamilyData;
use crate::App;

pub(super) fn simple_report(data: &FamilyData) {
    let output_pdf = "typst_report.pdf";
    let typst_template = work_on_me();

    // Here be dragons
    let typst_cache = App::DATA_DIR.to_string() + "/typst_cache";
    let pkg_resolver = PackageResolver::builder()
        .cache(FileSystemCache(PathBuf::from(typst_cache)))
        .build();

    let template = TypstEngine::builder()
        .main_file(typst_template)
        .search_fonts_with(TypstKitFontOptions::default())
        .add_file_resolver(pkg_resolver)
        .with_file_system_resolver(App::DATA_DIR)
        .build();

    let doc = template.compile().output
        .expect("typst::compile() returned an error!");

    let options = Default::default();

    let pdf = typst_pdf::pdf(&doc, &options).expect("Could not generate pdf.");
    std::fs::write(output_pdf, pdf).expect("Could not write pdf.");
}

fn work_on_me() -> String {
    let typst_template: &str = include_str!("./../assets/templates/report_result.typ");
    return typst_template.to_string();
}
