use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, checkbox, column, container, horizontal_rule, horizontal_space, pick_list, row, scrollable, text, text_input, tooltip};
use iced::widget::{Row, Column};
use iced::{Element, Length, Size, Padding};

use std::iter::zip;
use std::path::{Path, PathBuf};
use std::process::Command;
use native_dialog::FileDialog;
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
use std::collections::HashSet;
use std::error::Error;

use crate::{App, ContentPage, MotifFile};

#[derive(Debug, Clone)]
pub(crate) enum Message {
    Back,
    SetMotifs(MotifFile),
    MotifCheckbox(usize, bool),
    MotifGroupbox(usize, bool),

    BamChanged(String),
    SelectBam,
    EditMetadata(PathBuf, PathBuf),

    RunDante,
    OpenResults,
    Print,
    CheckboxOutBAM(bool),
}

#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
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

    motifs: Vec<(bool, String, Vec<String>, String)>,
    groups: Vec<(bool, String)>,
}

impl Data {
    pub(super) fn view(&self, size: Size) -> Element<Message> {
        let mut content = column![].align_x(Horizontal::Center);

        content = make_header(content);
        content = make_form(content, self);

        content = content.push(container(horizontal_rule(2)).padding(25));
        content = make_report(content, self, size);
        content = content.push(draw_open_button(self));

        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return scrollable(content).into();
    }

    pub(super) fn update(&mut self, m: Message) {
        match m {
            Message::BamChanged(content) => { self.bam_file = Some(PathBuf::from(content)); }
            Message::SelectBam => { load_file(&mut self.bam_file); },
            Message::RunDante => { run_analysis(self); },
            Message::OpenResults => { open_results(self); },
            Message::CheckboxOutBAM(is_checked) => self.out_bam = is_checked,

            Message::SetMotifs(motif_file) => { update_motif_selection(self, motif_file); },
            Message::MotifGroupbox(idx, checked) => { toggle_group(self, idx, checked); },
            Message::MotifCheckbox(idx, checked) => { self.motifs[idx].0 = checked; /* Task::none() */ },
            Message::Print => { println!(); }

            Message::Back => { unreachable!() },
            Message::EditMetadata(_, _) => { unreachable!() /* implemented in main */ }
        }
    }

    pub(super) fn init(path: PathBuf, analysis_name: String) -> ContentPage {
        let data = Data { path, analysis_name, ..Default::default() };
        data.save();
        ContentPage::AnalysisSingle(data)
    }

    fn save(&self) -> PathBuf {
        let json = serde_json::to_string(self).unwrap();
        let mut output = self.path.clone();
        output.push("params.json");
        let mut out = File::create(&output)
            .expect("Cannot open file for writing.");
        out.write_all(json.as_bytes())
            .expect("Cannot write to output file.");
        return output;
    }

    pub(super) fn load(mut path: PathBuf) -> Self {
        path.push("params.json");
        let json: String = std::fs::read_to_string(path)
            .expect("Cannot read file.");
        serde_json::from_str(&json)
            .expect("Cannot parse json.")
    }
}

fn toggle_group(data: &mut Data, idx: usize, checked: bool) {
    data.groups[idx].0 = checked;
    let group = data.groups[idx].1.clone();
    for x in &mut data.motifs { if x.2.contains(&group) { x.0 = checked; } }
}

fn update_motif_selection(data: &mut Data, motif_file: MotifFile) {
    match motif_file {
        MotifFile::Custom => {
            if let Ok(Some(path)) = FileDialog::new().show_open_single_file() {
                data.selected = Some(motif_file);
                data.selected_file = Some(path);

                let format = validate_STR_format(data.selected_file.as_ref().unwrap());
                if format.is_ok() {
                    data.motifs = parse_motifs(data.selected_file.as_ref().unwrap());
                    data.groups = get_groups(data.motifs.as_ref());
                    data.message_line = "".to_string();
                } else {
                    data.motifs = Vec::new();
                    data.groups = Vec::new();
                    data.message_line = format.unwrap_err().to_string();
                }
            }
        },
        _ => {
            let motif_str = motif_file.to_string();
            let path = PathBuf::from(App::DATA_DIR.to_string() + "/" + &motif_str + ".tsv");
            data.selected = Some(motif_file);
            data.selected_file = Some(path);
            data.motifs = parse_motifs(data.selected_file.as_ref().unwrap());
            data.groups = get_groups(data.motifs.as_ref());
            data.message_line = "".to_string();
        }
    }
}

#[allow(non_snake_case)]
fn validate_STR_format(path: &Path) -> Result<(), Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let header = reader.lines().next().ok_or("File does not contain any lines.")??;
    let header: Vec<_> = header.split('\t').collect();

    if header[0] != "Disease ID"
        { return Err("1st column has incorrect name".into()); }
    if header[1] != "HGVS nomenclature (GRCh38 reference)"
        { return Err("2nd column has incorrect name".into()); }
    if header[2] != "Left flank"
        { return Err("3rd column has incorrect name".into()); }
    if header[3] != "Right flank"
        { return Err("4th column has incorrect name".into()); }
    if header[4] != "Groups"
        { return Err("5th column has incorrect name".into()); }
    if header[5] != "Disease name"
        { return Err("6th column has incorrect name".into()); }

    return Ok(());
}

fn parse_motifs(path: &Path) -> Vec<(bool, String, Vec<String>, String)> {
    let file = File::open(path).expect("Cannot find motif file.");
    let reader = BufReader::new(file);

    let mut result = Vec::new();
    for line in reader.lines().skip(1) {
        let line = line.expect("Cannot read line from motif file.").trim().to_string();
        let split: Vec<_> = line.split('\t').collect();

        let id = split[0].to_string();
        let groups = split[4].split(',').map(|x| x.to_string()).collect();
        let description = split[5].to_string();
        result.push((false, id, groups, description));
    }

    return result;
}

fn get_groups(motifs: &[(bool, String, Vec<String>, String)]) -> Vec<(bool, String)> {
    let groups: HashSet<(bool, String)> = motifs.iter()
        .flat_map(|x| &x.2)
        .map(|x| (false, x.to_string()))
        .collect();
    let mut groups: Vec<_> = groups.into_iter().collect();
    groups.sort();
    return groups;
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

    let metadata = get_metadata(data.bam_file.clone());
    let edit_button = if data.bam_file.is_some() {
        let edit_msg = Message::EditMetadata(data.path.clone(), data.bam_file.clone().unwrap());
        button("Edit metadata").on_press(edit_msg)
    } else {
        button("Edit metadata")
    };

    row![
        container(text("BAM file: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(
            tooltip(
                text_input("Type path or click search...", &proband).on_input(Message::BamChanged),
                container(text(metadata)).padding(5).style(container::rounded_box),
                tooltip::Position::FollowCursor,
            )
        ).padding(App::PAD1),
        // container(pick_list(sex, data.proband_sex, Message::ProbandSetSex).placeholder("sex").width(PSIZE)).padding(App::PAD1),
        container(button("Search").on_press(Message::SelectBam)).padding(App::PAD1),
        container(edit_button).padding(App::PAD1)
    ].padding(10).align_y(Vertical::Center)
}

fn get_metadata(bam_file: Option<PathBuf>) -> String {
    if bam_file.is_none() { return "No metadata found.".to_string(); }

    let mut meta_file = bam_file.unwrap();
    meta_file.set_extension("meta.tsv");
    if !meta_file.exists() { return "No metadata found.".to_string(); }
    if meta_file.is_dir() { return "No metadata found.".to_string(); }

    let mut lines = BufReader::new(File::open(meta_file).expect("Cannot open metadata file.")).lines();
    let header = lines.next().unwrap().unwrap();
    let header = header.split("\t");
    let content = lines.next().unwrap().unwrap();
    let content = content.split("\t");

    // TODO: select most important data
    // let mut result = String::new();
    // for (h, c) in zip(header, content) {
    //     result.push_str(h);
    //     result.push(':');
    //     result.push_str(c);
    //     result.push('\n');
    // }
    let mut result = "".to_string();
    result.push_str("Metadata stored in vpuk-23-001504-A.meta.tsv\n");
    result.push_str("Patient name: John Doe\n");
    result.push_str("Sample ID: 18298371\n");
    result.push_str("Gender: Male\n");
    result.push_str("+ 32 other entries.\n");

    return result;
}

fn make_motif_selection(selected: Option<MotifFile>, selected_file: &Option<PathBuf>) -> Row<Message> {
    let motif_files = [MotifFile::STRSet_20250311, MotifFile::Custom];

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

fn run_analysis(data: &Data) {
    let out_dir = data.path.to_string_lossy().to_string();
    data.save();

    // required params
    let Some(ref bam_file) = data.bam_file else { return; };
    let Some(ref motif_file) = data.selected_file else { return; };
    let output: String = out_dir.clone() + "/remaSTR_result.tsv";

    // optional params
    let out_bam = data.out_bam;
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

fn make_report<'a>(mut content: Column<'a, Message>, data: &'a Data, size: Size) -> Column<'a, Message> {
    let available_width = size.width as usize - 5 - 160 - 5 - 5;

    let mut i = 0;
    let mut r = row![].padding(5).align_y(Vertical::Center);
    r = r.push(container(text("Group filter: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1));
    r = r.extend(make_group_row(&data.groups, available_width, &mut i));
    r = r.push(horizontal_space());
    // let r = std::convert::Into::<Element<Message>>::into(r).explain(iced::Color::BLACK);
    content = content.push(r);

    let mut i = 0;
    let mut r = row![].padding(5).align_y(Vertical::Center);
    r = r.push(container(text("Motif filter: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1));
    r = r.extend(make_checkbox_row(&data.motifs, available_width, &mut i));
    r = r.push(horizontal_space());
    // let r = std::convert::Into::<Element<Message>>::into(r).explain(iced::Color::BLACK);
    content = content.push(r);

    // BUG: if available_width is too small, i is never increased
    while i < data.motifs.len() {
        let mut r = row![].padding(5).align_y(Vertical::Center);
        r = r.push(container(text("")).width(160).align_x(Horizontal::Right));
        r = r.extend(make_checkbox_row(&data.motifs, available_width, &mut i));
        r = r.push(horizontal_space());
        // let r = std::convert::Into::<Element<Message>>::into(r).explain(iced::Color::BLACK);
        content = content.push(r);
    }

    const PAD2: Padding = Padding { bottom: 0.0, top: 10.0, right: 15.0, left: 0.0 };
    let r = row![
        container(text("")).width(160),
        container(button("View")).padding(PAD2),
        container(button("Print").on_press(Message::Print)).padding(PAD2),
        horizontal_space(),
    ].padding(10).align_y(Vertical::Center);
    // let r = std::convert::Into::<Element<Message>>::into(r).explain(iced::Color::BLACK);
    content = content.push(r);
    return content;
}

fn make_group_row<'a>(groups: &'a[(bool, String)], available_width: usize, i: &mut usize) -> Vec<Element<'a, Message>> {
    const PAD: Padding = Padding { bottom: 0.0, top: 0.0, right: 15.0, left: 0.0};
    let spacing = 15 /*checkbox*/ + 10 /*between checkbox and label*/ + 15 /*right padding*/;
    let letter_width = 11;

    let mut v = Vec::new();
    let mut cur_width = 0;
    while *i < (*groups).len() && cur_width + spacing + groups[*i].1.len() * letter_width < available_width {
        let (ref checked, ref id) = &groups[*i];
        let ii = *i;
        let f = move |b| Message::MotifGroupbox(ii, b);
        v.push(container(checkbox(id, *checked).on_toggle(f)).padding(PAD).into());
        cur_width += spacing + id.len() * letter_width;
        *i += 1;
    }
    return v;
}

fn make_checkbox_row<'a>(motifs: &'a[(bool, String, Vec<String>, String)], available_width: usize, i: &mut usize) -> Vec<Element<'a, Message>> {
    // TODO: are the lifetimes correct?
    const PAD: Padding = Padding { bottom: 0.0, top: 0.0, right: 15.0, left: 0.0 };
    let spacing = 15 /*checkbox*/ + 10 /*between checkbox and label*/ + 15 /*right padding*/;
    let letter_width = 11;

    let mut v = Vec::new();
    let mut cur_width = 0;
    while *i < (*motifs).len() && cur_width + spacing + motifs[*i].1.len() * letter_width < available_width {
        let (ref checked, ref id, _, ref name) = &motifs[*i];
        let ii = *i;
        let f = move |b| Message::MotifCheckbox(ii, b);
        v.push(container(
            tooltip(
                checkbox(id, *checked).on_toggle(f),
                container(text(name)).padding(5).style(container::rounded_box),
                tooltip::Position::FollowCursor,
            )
        ).padding(PAD).into());
        cur_width += spacing + id.len() * letter_width;
        *i += 1;
    }
    return v;
}


