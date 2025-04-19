use chrono::NaiveDate;
use icu_locid::locale;
use spreadsheet_ods::color::Rgb;
use spreadsheet_ods::formula;
use spreadsheet_ods::inch;
use spreadsheet_ods::mm;
use spreadsheet_ods::style::units::Margin;
use spreadsheet_ods::style::CellStyle;
use spreadsheet_ods::style::units::{Border, TextRelief};
use spreadsheet_ods::format;
use spreadsheet_ods::style::MasterPage;
use spreadsheet_ods::style::PageStyle;
use spreadsheet_ods::Length;
use spreadsheet_ods::{Sheet, Value, WorkBook};
use std::path::Path;

fn create_empty_workbook() -> WorkBook {
    let mut wb = WorkBook::new(locale!("en_US"));
    let mut sheet = Sheet::new("Dante - One-page report");
    sheet.set_value(0, 0, true);
    wb.push_sheet(sheet);
    let mut sheet = Sheet::new("Dante - Summary report");
    sheet.set_value(0, 0, true);
    wb.push_sheet(sheet);
    let mut sheet = Sheet::new("Dante - Results report");
    sheet.set_value(0, 0, true);
    wb.push_sheet(sheet);
    let mut sheet = Sheet::new("Dante - Technical report");
    sheet.set_value(0, 0, true);
    wb.push_sheet(sheet);
    return wb;
}

fn get_sheet_height_in_mm(sheet: &Sheet, max: u32) -> f64 {
    let mut result = 0.0;
    for i in 0..max {
        let x = sheet.row_height(i);
        print!("{:?} ", x);
        match x {
            Length::Default => { println!("Contains default length.") },
            Length::Cm(x) => { result += 10.0 * x; },
            Length::Mm(x) => { result += x; },
            Length::In(x) => { result += 25.4 * x; /* 1inch is exactly 25.4mm */ },
            Length::Pt(x) => { result += x / 72.0 * 25.4 /* 1pt is 1/72 of inch */ },
            Length::Pc(x) => { result += x * 12.0 / 72.0 * 25.4 /* 1pc is 12pt */},
            Length::Em(_) => { println!("Contains em length.") },
        }
    }
    println!();
    return result;
}

fn get_sheet_width_in_mm(sheet: &Sheet, max: u32) -> f64 {
    let mut result = 0.0;
    for i in 0..max {
        let x = sheet.col_width(i);
        print!("{:?} ", x);
        match x {
            Length::Default => { println!("Contains default length.") },
            Length::Cm(x) => { result += 10.0 * x; },
            Length::Mm(x) => { result += x; },
            Length::In(x) => { result += 25.4 * x; /* 1inch is exactly 25.4mm */ },
            Length::Pt(x) => { result += x / 72.0 * 25.4 /* 1pt is 1/72 of inch */ },
            Length::Pc(x) => { result += x * 12.0 / 72.0 * 25.4 /* 1pc is 12pt */},
            Length::Em(_) => { println!("Contains em length.") },
        }
    }
    println!();
    return result;
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
        for i in 1..=24 {
            sheet.clear_colstyle(i);
            sheet.set_col_width(i, mm!(8));
        }
        sheet.clear_colstyle(25);
        sheet.set_col_width(25, mm!(6));
    }
}

fn set_print_styles(wb: &mut WorkBook) {
    let mut ps = PageStyle::new_empty();
    ps.set_name("Printing");
    ps.set_page_width(mm!(210));
    ps.set_page_height(mm!(297));
    ps.set_margin_top(mm!(4.5));
    ps.set_margin_bottom(mm!(4.5));
    ps.set_margin_left(mm!(3));
    ps.set_margin_right(mm!(0));
    // ps.set_table_centering(center);
    ps.headerstyle_mut().set_height(mm!(0));
    ps.footerstyle_mut().set_height(mm!(0));

    let mut mpage: MasterPage = MasterPage::new_empty();
    mpage.set_name("Default".to_string());
    mpage.set_pagestyle(&ps.style_ref());
    wb.add_masterpage(mpage);

    wb.add_pagestyle(ps);
}

fn main() {
    let inpath = Path::new("template.ods");
    let outpath = Path::new("template2.ods");
    println!("A4 is 210 x 297 mm.");

    let mut wb = if inpath.exists() {
        println!("Reading existing file {:?}.", inpath);
        spreadsheet_ods::read_ods(inpath).unwrap()
    } else {
        println!("Creating new file.");
        create_empty_workbook()
    };

    let n = wb.num_sheets();
    println!("Template has {n} sheets.");

    for i in 0..n {
        let sheet = wb.sheet(i);
        let (n_rows, n_cols) = sheet.used_grid_size();

        let height = get_sheet_height_in_mm(sheet, n_rows);
        let width = get_sheet_width_in_mm(sheet, n_cols);

        println!("{}. {}", i, sheet.name());
        println!("\trows: {}, cols: {}", n_rows, n_cols);
        println!("\t{width:.3} x {height:.3} mm") // A4 is 210 x 297 mm
    }

    let x = wb.iter_masterpages().count();
    println!("# of masterpages: {x}");
    for x in wb.iter_masterpages() {
        println!("{}", x.name());
    }

    let x = wb.iter_pagestyles().count();
    println!("# of pagestyles: {x}");
    for x in wb.iter_pagestyles() {
        println!("{}", x.name());
    }

    normalize_rows_and_cols(&mut wb);
    set_print_styles(&mut wb);
    spreadsheet_ods::write_ods(&mut wb, outpath).expect("Failed to write ODS file.")
}

// fn fn1() {
//     let page_layout = PageLayoutBuilder::new("CustomLayout")
//         .margin_left(Length::cm(2.0))
//         .margin_right(Length::cm(2.0))
//         .margin_top(Length::cm(2.5))
//         .margin_bottom(Length::cm(2.5))
//         .build();
// 
//     // Assign layout to a page master
//     let page_master = PageMaster::new("pm1", "CustomLayout");
// 
//     // Add the page style using the page master
//     let page_style = PageStyle::new("MyPageStyle", "pm1");
// 
//     let mut book = WorkBook::new("WithMargins");
//     let mut sheet = Sheet::new("Sheet1");
// 
//     sheet.set_value((0, 0), "Hello, margins!");
// 
//     book.add_page_layout(page_layout);
//     book.add_page_master(page_master);
//     book.add_page_style(page_style);
// 
//     // Assign page style to the sheet
//     sheet.set_page_style("MyPageStyle");
// 
//     book.push_sheet(sheet);
// 
//     book.save("margins.ods").expect("Failed to save");
// }

        // for (pos, c) in sheet.iter_cols((0, 0)..(grid.0, grid.1)) {
        //     // println!("{:?}", c);
        //     // println!("{:?}", pos);
        //     // println!("{:?}", c.matrix_row_span());
        // }


    // if wb.num_sheets() == 0 {
    //     let mut sheet = Sheet::new("one");
    //     sheet.set_value(0, 0, true);
    //     wb.push_sheet(sheet);
    // }

    // let sheet = wb.sheet(0);
    // let _n = sheet.value(0, 0).as_f64_or(0f64);
    // if let Value::Boolean(v) = sheet.value(1, 1) {
    //     if *v { println!("was true"); }
    // }

    // if wb.num_sheets() == 1 { wb.push_sheet(Sheet::new("two")); }

    // let date_format = format::create_date_dmy_format("date_format");
    // let date_format = wb.add_datetime_format(date_format);

    // let mut date_style = CellStyle::new("nice_date_style", &date_format);
    // date_style.set_font_bold();
    // date_style.set_font_relief(TextRelief::Engraved);
    // date_style.set_border(mm!(0.2), Border::Dashed, Rgb::new(192, 72, 72));
    // let date_style_ref = wb.add_cellstyle(date_style);

    // let sheet = wb.sheet_mut(1);
    // sheet.set_value(0, 0, 21.4f32);
    // sheet.set_value(0, 1, "foo");
    // sheet.set_styled_value(0, 2, NaiveDate::from_ymd_opt(2020, 3, 1), &date_style_ref);
    // sheet.set_formula(0, 3, format!("of:={}+1", formula::fcellref(0, 0)));

    // let mut sheet = Sheet::new("sample");
    // sheet.set_value(5, 5, "sample");
    // wb.push_sheet(sheet);


