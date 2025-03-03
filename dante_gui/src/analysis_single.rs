use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, checkbox, column, container, horizontal_rule, horizontal_space, pick_list, row, scrollable, text, text_input};
use iced::widget::{Row, Column};
use iced::{Element, Size, Length};
use std::path::PathBuf;
use std::env;
use std::path::Path;
use std::process::Command;
use remastr;
use native_dialog::FileDialog;
use std::fs;

use crate::{App, ContentPage, MotifFile};

#[derive(Debug, Clone)]
pub(crate) enum Message {
    Back,
    SetMotifs(MotifFile),
    // MotifCheckbox(usize, bool),
    // MotifGroupbox(usize, bool),

    BamChanged(String),
    SelectBam,

    RunDante,
    OpenResults,
    CheckboxOutBAM(bool),
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct Data {
    path: PathBuf,
    analysis_name: String,
    selected: Option<MotifFile>,
    selected_file: Option<PathBuf>,

    bam_file: Option<PathBuf>,
    motif_file: Option<PathBuf>,
    output: Option<PathBuf>,
    out_bam: bool,
    message_line: String,
}

impl Data {
    pub(super) fn init(path: PathBuf, analysis_name: String) -> ContentPage {
        ContentPage::AnalysisSingle(Data {
            path, analysis_name, ..Default::default()
        })
    }

    pub(super) fn view(&self, size: Size) -> Element<Message> {
        println!("{:?}", size);
        let mut content = column![].align_x(Horizontal::Center);

        content = make_header(content);
        content = make_form(content, self);

        content = content.push(container(horizontal_rule(2)).padding(25));
        // content = make_report(content, self, size);
        content = content.push(draw_open_button(self));

        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return scrollable(content).into();
    }

    pub(super) fn update(&mut self, m: Message) {
        match m {
            Message::BamChanged(content) => { self.bam_file = Some(PathBuf::from(content)); }
            Message::SelectBam => { load_file(&mut self.bam_file); },
            Message::RunDante => { run1(self); },
            Message::OpenResults => { open_results(self); },
            Message::CheckboxOutBAM(is_checked) => self.out_bam = is_checked,

            Message::SetMotifs(motif_file) => { update_motif_selection(self, motif_file); },
            // Message::MotifGroupbox(idx, name) => { println!("{} {}", idx, name); },
            // Message::MotifCheckbox(idx, name) => { println!("{} {}", idx, name); },

            Message::Back => { unreachable!() },
        }
    }
}

fn update_motif_selection(data: &mut Data, motif_file: MotifFile) {
    match motif_file {
        MotifFile::STRSet_20220902 => {
            let path = PathBuf::from(App::DATA_DIR.to_string() + "/STRSet_20220902.tsv");
            data.selected = Some(motif_file);
            data.selected_file = Some(path);
            // data.motifs = parse_motifs(data.selected_file.as_ref().unwrap());
            // data.groups = get_groups(data.motifs.as_ref());
        },
        MotifFile::Custom => {
            if let Ok(Some(path)) = FileDialog::new().show_open_single_file() {
                data.selected = Some(motif_file);
                data.selected_file = Some(path);
                // TODO: How to handle incorrect inputs?
                // data.motifs = parse_motifs(data.selected_file.as_ref().unwrap());
                // data.groups = get_groups(data.motifs.as_ref());
            }
        },
    }
}

fn make_header(mut content: Column<Message>) -> Column<Message> {
    content = content.push(row![
        container(button("Back").on_press(Message::Back)).width(100),
        container(text("Single analysis").size(App::H1_SIZE)).align_x(Horizontal::Center).width(Length::Fill),
        container("").width(100),
    ].padding(25).align_y(Vertical::Center));
    return content;
}

fn make_form<'a>(mut content: Column<'a, Message>, data: &'a Data) -> Column<'a, Message> {
    content = content.push(row![
        container(text("Analysis name: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text(data.analysis_name.clone())).width(Length::Fill).align_x(Horizontal::Left)
    ].padding(10).align_y(Vertical::Center));

    content = content.push(make_motif_selection(data.selected, &data.selected_file));
    content = content.push(make_proband_row(data));

    content = content.push(row![
        container(text("")).width(160),
        container(checkbox("Output filtered BAM", data.out_bam).on_toggle(Message::CheckboxOutBAM)),
        horizontal_space(),
    ].padding(10.0).align_y(Vertical::Center));

    content = content.push(row![
        container("").width(160),
        button("Run").on_press(Message::RunDante),
        container(text(data.message_line.clone()).align_x(Horizontal::Left)).padding(App::PAD2),
        horizontal_space(),
    ].padding(10.0).align_y(Vertical::Center));

    return content;
}

fn make_proband_row(data: &Data) -> Row<Message> {
    let proband = data.bam_file.clone().unwrap_or_default().to_string_lossy().to_string();
    // let sex = [Sex::Male, Sex::Female, Sex::Unknown];

    row![
        container(text("BAM file: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text_input("Type path or click search...", &proband).on_input(Message::BamChanged)).padding(App::PAD1),
        // container(pick_list(sex, data.proband_sex, Message::ProbandSetSex).placeholder("sex").width(PSIZE)).padding(App::PAD1),
        container(button("Search").on_press(Message::SelectBam)).padding(App::PAD1)
    ].padding(10).align_y(Vertical::Center)
}

fn make_motif_selection(selected: Option<MotifFile>, selected_file: &Option<PathBuf>) -> Row<Message> {
    let motif_files = [MotifFile::STRSet_20220902, MotifFile::Custom];

    let content = match selected {
        None => row![
            container(text("Motifs: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
            container(pick_list(motif_files, selected, Message::SetMotifs).placeholder("type")),
            container(text("")).width(Length::Fill).align_x(Horizontal::Left)
        ].padding(10).align_y(Vertical::Center),
        Some(MotifFile::Custom) => {
            let x: String = selected_file.clone().unwrap().to_string_lossy().to_string();
            row![
                container(text("Motifs: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
                container(pick_list(motif_files, selected, Message::SetMotifs).placeholder("type")).padding(App::PAD1),
                container(text(x)).width(Length::Fill).align_x(Horizontal::Left)
            ].padding(10).align_y(Vertical::Center)
        }
        Some(_) => {
            row![
                container(text("Motifs: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
                container(pick_list(motif_files, selected, Message::SetMotifs).placeholder("type")),
                container(text("")).width(Length::Fill).align_x(Horizontal::Left)
            ].padding(10).align_y(Vertical::Center)
        }
    };
    return content;
}

fn run1(state: &Data) {
    println!("{:?}", state);

    // required params
    let Some(ref bam_file) = state.bam_file else { return; };
    let Some(ref motif_file) = state.selected_file else { return; };
    let out_dir = state.path.to_string_lossy().to_string();
    let output: String = out_dir.clone() + "/remaSTR_result.tsv";

    // optional params
    let out_bam = state.out_bam;
    let dedup = false;
    let print_quality = false;
    let q = 30;
    let score: Option<char> = None;

    remastr::run(bam_file, motif_file, output.clone(), out_bam, (dedup, q, score, print_quality));
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
    let output = state.path.to_string_lossy().to_string() + "/report.html";
    opener::open(output).unwrap();
}

fn draw_open_button<'a>(state: &Data) -> Element<'a, Message> {
    let output = state.path.to_string_lossy().to_string() + "/report.html";
    if Path::new(&output).exists() { 
        let report_line = format!("Report file stored in {}.", output);
        row![
            container("").width(160),
            button("Open results").on_press(Message::OpenResults),
            container(text(report_line).align_x(Horizontal::Left)).padding(App::PAD2),
            horizontal_space(),
        ].padding(10.0).align_y(Vertical::Center).into()
    } else {
        let report_line = "No report file present.".to_string();
        row![
            container("").width(160),
            button("Open results"),
            container(text(report_line).align_x(Horizontal::Left)).padding(App::PAD2),
            horizontal_space(),
        ].padding(10.0).align_y(Vertical::Center).into()
    }
}

fn load_file(result: &mut Option<PathBuf>) {
    let path = FileDialog::new().show_open_single_file().unwrap();
    let path = match path {
        Some(path) => path,
        None => return,
    };
    *result = Some(path);
}
