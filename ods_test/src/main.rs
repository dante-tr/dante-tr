use chrono::NaiveDate;
use icu_locid::locale;
use spreadsheet_ods::color::Rgb;
use spreadsheet_ods::formula;
use spreadsheet_ods::inch;
use spreadsheet_ods::mm;
use spreadsheet_ods::style::CellStyle;
use spreadsheet_ods::style::units::{Border, TextRelief};
use spreadsheet_ods::format;
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

fn main() {
    let inpath = Path::new("template.ods");
    let outpath = Path::new("template.ods");

    let mut wb = if inpath.exists() {
        println!("Reading existing file {:?}.", inpath);
        spreadsheet_ods::read_ods(inpath).unwrap()
    } else {
        println!("Creating new file.");
        create_empty_workbook()
    };

    let n = wb.num_sheets();
    println!("Template has {n} sheets.\n");

    for i in 0..n {
        let sheet = wb.sheet_mut(i);
        print!("{}. {} - ", i, sheet.name());
        let grid = sheet.used_grid_size();
        println!("rows: {}, cols: {}", grid.0, grid.1);
        let mut width = 0.0;
        for i in 0..grid.1 {
            let x = sheet.col_width(i);
            let Length::In(x) = x else { todo!() };
            width += x;
        }
        println!("Width = {width} inches.");

        let mut height = 0.0;
        for i in 0..grid.0 {
            sheet.set_row_height(i, inch!(1.0));
            let x = sheet.row_height(i);
            println!("{i} set to {:?}", x);
            // let Length::In(x) = x else { todo!() };
            // height += x;
        }
        println!("Height = {height} inches.");
        // A4 is 210 x 297 mm
        // 1 inch is exactly 25.4mm
    }
    spreadsheet_ods::write_ods(&mut wb, outpath).expect("Failed to write ODS file.")
}

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


