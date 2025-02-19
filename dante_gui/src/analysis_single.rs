use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, checkbox, column, container, horizontal_rule, row, text, text_input};
use iced::Element;
use std::path::PathBuf;
use std::env;
use std::path::Path;
use std::process::Command;
use remastr::run;
use native_dialog::FileDialog;
use std::fs;

use crate::App;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Data {
    bam_file: Option<PathBuf>,
    motif_file: Option<PathBuf>,
    output: Option<PathBuf>,
    out_bam: bool,
    message_line: String,
}

#[derive(Debug, Clone)]
pub(crate) enum Message {
    Back,
    BamChanged(String),
    SelectBam,
    MotifChanged(String),
    SelectMotif,
    OutdirChanged(String),
    SelectOutdir,
    RunDante,
    OpenResults,
    CheckboxOutBAM(bool),
}

pub(crate) fn update(data: &mut Data, m: Message) {
    match m {
        Message::BamChanged(content) => { data.bam_file = Some(PathBuf::from(content)); }
        Message::MotifChanged(content) => { data.motif_file = Some(PathBuf::from(content)); },
        Message::OutdirChanged(content) => { data.output = Some(PathBuf::from(content)); },
        Message::SelectBam => { load_file(&mut data.bam_file); },
        Message::SelectMotif => { load_file(&mut data.motif_file); },
        Message::SelectOutdir => { load_dir(&mut data.output); },
        Message::RunDante => { run1(data); },
        Message::OpenResults => { open_results(data); },
        Message::CheckboxOutBAM(is_checked) => data.out_bam = is_checked,
        Message::Back => { unreachable!() },
    }
}

pub fn view<'a>(state: &'a App, data: &'a Data) -> Element<'a, Message> {
    column![
        button("Back").on_press(Message::Back),
        loader_row("BAM file:", &data.bam_file, Message::BamChanged, Message::SelectBam),
        loader_row("Motif file:", &data.motif_file, Message::MotifChanged, Message::SelectMotif),
        horizontal_rule(2),
        loader_row("Output directory:", &data.output, Message::OutdirChanged, Message::SelectOutdir),

        row![
            container("").width(App::LEFT_WIDTH).padding(App::PAD1),
            checkbox("Output BAM", data.out_bam).on_toggle(Message::CheckboxOutBAM),
        ].padding(10.0).align_y(Vertical::Center),
        row![
            container("").width(App::LEFT_WIDTH).padding(App::PAD1),
            button("Run").on_press(Message::RunDante),
            container(text(data.message_line.clone()).align_x(Horizontal::Left)).padding(App::PAD2),
        ].padding(10.0).align_y(Vertical::Center),
        draw_open_button(data),
    ].align_x(Horizontal::Left).into()
}

fn draw_open_button<'a>(state: &Data) -> Element<'a, Message> {
    let report_present;
    let report_line;
    match &state.output {
        Some(x) => {
            let mut x = path_to_string(x);
            x.push_str("/report.html");
            if Path::new(&x).exists() {
                report_present = true;
                report_line = format!("Report file stored in {}.", x);
            } else {
                report_present = false;
                report_line = "No report file present.".to_string();
            }
        },
        None => {
            report_present = false;
            report_line = "No report file present.".to_string();
        }
    };

    if report_present {
        row![
            container("").width(App::LEFT_WIDTH).padding(App::PAD1),
            button("Open results").on_press(Message::OpenResults),
            container(text(report_line).align_x(Horizontal::Left)).padding(App::PAD2),
        ].padding(10.0).align_y(Vertical::Center).into()
    } else {
        row![
            container("").width(App::LEFT_WIDTH).padding(App::PAD1),
            button("Open results"),
            container(text(report_line).align_x(Horizontal::Left)).padding(App::PAD2),
        ].padding(10.0).align_y(Vertical::Center).into()
    }
}

fn loader_row<'a>(
    desc: &'a str, filename: &'a Option<PathBuf>, on_input: impl Fn(String) -> Message + 'a, on_press: Message
) -> Element<'a, Message> {
    let filename_str: String = match filename.as_ref() {
        Some(x) => path_to_string(x),
        None => "".to_string()
    };

    row![
        container(text(desc).width(App::LEFT_WIDTH).align_x(Horizontal::Right)).padding(App::PAD1),
        text_input("Type path or click search...", &filename_str).on_input(on_input).font(App::BOLD_MONO),
        container(button("Search").on_press(on_press)).padding(App::PAD2),
    ].padding(10.0).align_y(Vertical::Center).into()
}

fn path_to_string(path: &Path) -> String {
    let cwd = env::current_dir().unwrap().display().to_string();
    match path.strip_prefix(cwd) {
        Ok(x) => { x.display().to_string() },
        Err(_) => { path.display().to_string() }
    }
}

fn run1(state: &Data) {
    println!("{:?}", state);

    // required params
    let Some(ref bam_file) = state.bam_file else {
        return;
    };
    let Some(ref motif_file) = state.motif_file else {
        return;
    };

    let mut output: String = match state.output {
        Some(ref x) => x.display().to_string(),
        None => {
            return;
        },
    };
    let out_dir = output.clone();
    if !Path::new(&out_dir).exists() {
        fs::create_dir(&out_dir).expect("Cannot create directory.");
    }
    output.push_str("/remaSTR_result.tsv");

    // optional params
    let out_bam = state.out_bam;
    let dedup = false;
    let print_quality = false;
    let q = 30;
    let score: Option<char> = None;

    run(bam_file, motif_file, output.clone(), out_bam, (dedup, q, score, print_quality));
    println!("remaSTR finished.");
    // self.message_line = "remaSTR finished.".to_string();
    let bin = format!("{}/dante_remastr_standalone", App::DATA_DIR);
    let output_log = Command::new(bin)
        .arg("--input-tsv").arg(output.clone())
        .arg("--output-dir").arg(out_dir)
        .arg("--verbose")
        .output()
        .expect("failed to run python part of Dante");
    println!("Dante finished.");
    // self.message_line = "Dante finished.".to_string();
    println!("{:?}", output_log);
}

fn open_results(state: &mut Data) {
    match state.output.as_ref() {
        Some(x) => {
            let mut output: String = x.to_str().unwrap().to_string();
            output.push_str("/report.html");
            opener::open(output).unwrap();
        },
        None => {
            state.message_line = "No report found.".to_string();
        },
    }
}

fn load_file(result: &mut Option<PathBuf>) {
    // TODO: "." does not work under Windows
    // let path = FileDialog::new().set_location(".").show_open_single_file().unwrap();
    let path = FileDialog::new().show_open_single_file().unwrap();
    let path = match path {
        Some(path) => path,
        None => return,
    };
    *result = Some(path);
}

fn load_dir(result: &mut Option<PathBuf>) {
    // TODO: "." does not work under Windows
    // let path = FileDialog::new().set_location(".").show_open_single_dir().unwrap();
    let path = FileDialog::new().show_open_single_dir().unwrap();
    let path = match path {
        Some(path) => path,
        None => return,
    };
    *result = Some(path);
}
