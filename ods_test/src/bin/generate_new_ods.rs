use icu_locid::locale;
use spreadsheet_ods::style::units::PrintCentering;
use spreadsheet_ods::{mm, pt, Value};
use spreadsheet_ods::style::MasterPage;
use spreadsheet_ods::style::PageStyle;
use spreadsheet_ods::{Sheet, WorkBook};
use std::path::Path;
use std::env;

const WIDTH: u32 = 24;
const ROWS_PER_PAGE: u32 = 48;

fn create_new_workbook() -> WorkBook {
    let mut wb = WorkBook::new(locale!("en_US"));

    let mut sheet = Sheet::new("Dante - Results report");
    for i in 0..WIDTH { sheet.set_value(0, i, i); }
    for j in 0..ROWS_PER_PAGE { sheet.set_value(j, 0, j); }
 
    // let value: Value;
    // sheet.set_value(1, 1, value);
    wb.push_sheet(sheet);

    let mut sheet = Sheet::new("Dante - Summary report");
    for i in 0..WIDTH { sheet.set_value(0, i, i); }
    for j in 0..ROWS_PER_PAGE { sheet.set_value(j, 0, j); }
    wb.push_sheet(sheet);

    let mut sheet = Sheet::new("Dante - One-page report");
    for i in 0..WIDTH { sheet.set_value(0, i, i); }
    for j in 0..ROWS_PER_PAGE { sheet.set_value(j, 0, j); }
    wb.push_sheet(sheet);

    // let mut sheet = Sheet::new("Dante - Technical report");
    // for i in 0..WIDTH { sheet.set_value(0, i, i); }
    // for j in 0..ROWS_PER_PAGE { sheet.set_value(j, 0, j); }
    // wb.push_sheet(sheet);

    return wb;
}

fn normalize_rows_and_cols(wb: &mut WorkBook) {
    let n = wb.num_sheets();
    for k in 0..n {
        let sheet = wb.sheet_mut(k);
        let (n_rows, _) = sheet.used_grid_size();

        for i in 0..n_rows {
            sheet.clear_rowstyle(i);
            sheet.set_row_height(i, mm!(6));
        }

        sheet.clear_colstyle(0);
        sheet.set_col_width(0, mm!(6));
        for i in 1..WIDTH-1 {
            sheet.clear_colstyle(i);
            // sheet.set_col_width(i, pt!(16));
            sheet.set_col_width(i, mm!(8));
        }
        sheet.clear_colstyle(WIDTH-1);
        sheet.set_col_width(WIDTH-1, mm!(6));
    }
}

fn set_print_styles(wb: &mut WorkBook) {
    let mut ps = PageStyle::new_empty();
    ps.set_name("Printing");
    ps.set_page_width(mm!(210));
    ps.set_page_height(mm!(297));
    ps.set_margin_top(pt!(21));
    ps.set_margin_bottom(pt!(21));
    ps.set_table_centering(PrintCentering::Horizontal);
    ps.set_margin_left(pt!(0));
    ps.set_margin_right(pt!(0));
    ps.headerstyle_mut().set_height(pt!(0));
    ps.footerstyle_mut().set_height(pt!(0));

    let mut mpage: MasterPage = MasterPage::new_empty();
    mpage.set_name("Default".to_string());
    mpage.set_pagestyle(&ps.style_ref());
    wb.add_masterpage(mpage);

    wb.add_pagestyle(ps);
}

fn main() {
    let outpath = env::args().nth(1).unwrap();
    let outpath = Path::new(&outpath);
    let mut wb = create_new_workbook();

    normalize_rows_and_cols(&mut wb);
    set_print_styles(&mut wb);
    spreadsheet_ods::write_ods(&mut wb, outpath).expect("Failed to write ODS file.")
}
