use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, checkbox, column, container, horizontal_rule, horizontal_space, pick_list, row, scrollable, text, text_input, tooltip};
use iced::widget::{Row, Column, Button};
use iced::{Element, Length, Size, Padding, Task};

use std::iter::zip;
use std::path::{Path, PathBuf};
use native_dialog::FileDialog;
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
use std::error::Error;

use crate::{App, ContentPage, MotifFile};
use crate::async_tasks;

#[derive(Debug, Clone)]
pub(crate) enum Message {
    Back,
    SetMotifs(MotifFile),
    MotifCheckbox(usize, bool),
    MotifGroupbox(usize, bool),

    BamChanged(String),
    SelectBam,
    EditMetadata(PathBuf, PathBuf),
    EditResults(Data),

    RunDante,
    AnalysisProgress(String),
    OpenResults,
    Print,
    SetReport(ReportType),
    CheckboxOutBAM(bool),
}

#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub(super) struct Data {
    analysis_path: PathBuf,
    analysis_name: String,

    selected: Option<MotifFile>,
    selected_file: Option<PathBuf>,

    bam_file: Option<PathBuf>,
    motif_file: Option<PathBuf>,
    output: Option<PathBuf>,
    out_bam: bool,
    message_line: String,

    selected_report: Option<ReportType>,
    message_line2: String,

    motifs: Vec<(bool, String, Vec<String>, String)>,  // (checked, id, groups, description)
    groups: Vec<(bool, String)>,
}

impl Data {
    pub(super) fn view(&self, size: Size) -> Element<Message> {
        let mut content = column![].align_x(Horizontal::Center);

        content = view_header(content);
        content = view_form(content, self);

        content = content.push(container(horizontal_rule(2)).padding(25));
        content = view_report(content, self, size);
        // content = content.push(draw_open_button(self));

        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return scrollable(content).into();
    }

    pub(super) fn update(&mut self, m: Message) -> Task<Message> {
        match m {
            Message::BamChanged(content) => { self.bam_file = Some(PathBuf::from(content)); Task::none() }
            Message::SelectBam => { load_file(&mut self.bam_file); Task::none() },
            Message::RunDante => { analyze(self) },
            Message::OpenResults => { open_results(self); Task::none() },
            Message::CheckboxOutBAM(is_checked) => { self.out_bam = is_checked; Task::none() },

            Message::SetMotifs(motif_file) => { update_motif_selection(self, motif_file); Task::none() },
            Message::MotifGroupbox(idx, checked) => { toggle_group(self, idx, checked); Task::none() },
            Message::MotifCheckbox(idx, checked) => { self.motifs[idx].0 = checked; Task::none() },
            Message::AnalysisProgress(msg) => { self.message_line = msg; Task::none() }
            Message::SetReport(report) => { self.selected_report = Some(report); Task::none() }
            Message::Print => { print_report(self); Task::none() }

            Message::Back => { unreachable!("Implemented in App::update."); },
            Message::EditMetadata(_, _) => { unreachable!("Implemented in App::update."); }
            Message::EditResults(_) => { unreachable!("Implemented in App::update.") }
        }
    }

    pub(super) fn init(path: PathBuf, analysis_name: String) -> ContentPage {
        let data = Data { analysis_path: path, analysis_name, ..Default::default() };
        data.save();
        ContentPage::AnalysisSingle(data)
    }

    pub(super) fn save(&self) -> PathBuf {
        let json = serde_json::to_string(self).unwrap();
        let mut output = self.analysis_path.clone();
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

    pub(super) fn get_checked_motif_ids(&self) -> Vec<String> {
        return self.motifs.iter().filter(|x| x.0).map(|x| x.1.replace("/", "_")).collect();
    }

    pub(crate) fn get_source(&self) -> PathBuf {
        return self.analysis_path.clone();
    }

    pub(crate) fn get_sample(&self) -> String {
        let Some(ref filepath) = self.bam_file else { panic!("There is no BAM file."); };
        let mut samplepath = filepath.clone();
        samplepath.set_extension("");
        return samplepath.file_name().expect("No filename.").to_string_lossy().to_string();
    }
}

fn toggle_group(data: &mut Data, idx: usize, checked: bool) {
    data.groups[idx].0 = checked;
    let group = data.groups[idx].1.clone();
    for x in &mut data.motifs { if x.2.contains(&group) { x.0 = checked; } }
}

fn update_motif_selection(data: &mut Data, motif_file: MotifFile) {
    use crate::analysis_common::{parse_motifs, get_groups};
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

fn view_header(mut content: Column<Message>) -> Column<Message> {
    content = content.push(row![
        container(button("Back").on_press(Message::Back)).width(100),
        container(text("Single analysis").size(App::H1_SIZE)).align_x(Horizontal::Center).width(Length::Fill),
        container("").width(100),
    ].padding(25).align_y(Vertical::Center));
    return content;
}

fn view_form<'a>(mut content: Column<'a, Message>, data: &'a Data) -> Column<'a, Message> {
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
    let (metadata, edit_button) = get_metadata(data);

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

fn get_metadata(data: &Data) -> (String, Button<Message>) {
    let bam_file = match &data.bam_file {
        None => { return ("No BAM file found.".to_string(), button("Edit metadata")); },
        Some(bam_file) => { bam_file },
    };
    if !bam_file.exists() { return ("No BAM file found.".to_string(), button("Edit metadata")); }

    let mut meta_file = bam_file.clone();
    meta_file.set_extension("meta.tsv");
    let edit_msg = Message::EditMetadata(data.analysis_path.clone(), meta_file.clone());
    if !meta_file.exists() {
        return ("No metadata found.".to_string(), button("Edit metadata").on_press(edit_msg));
    }

    let mut lines = BufReader::new(File::open(&meta_file).expect("Cannot open metadata file.")).lines();
    let header = lines.next().unwrap().unwrap();
    let header = header.split("\t");
    let content = lines.next().unwrap().unwrap();
    let content = content.split("\t");

    let mut metadata = String::new();
    let mut n_others = 0;
    metadata.push_str(&format!("Metadata stored in {}\n", meta_file.file_name().unwrap().to_str().unwrap()));
    for (h, c) in zip(header, content) {
        if !c.is_empty() {
            match h.strip_prefix("*") {
                Some(stripped) => {metadata.push_str(&format!("{}: {}\n", stripped, c)); },
                None => { n_others += 1; }
            }
        }
    }
    metadata.push_str(&format!("+ {} other entries.\n", n_others));

    return (metadata, button("Edit metadata").on_press(edit_msg));
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

fn analyze(data: &mut Data) -> Task<Message> {
    data.save();

    let Some(ref motif_file) = data.selected_file else { unreachable!() }; 
    let Some(ref bam_file) = data.bam_file else { unreachable!() };

    let mut output_file = data.analysis_path.clone();
    output_file.push(bam_file.file_name().unwrap());
    output_file.set_extension("");
    let dante_output_dir = output_file.clone();
    output_file.push("annotations.tsv");

    if !dante_output_dir.exists() {
        std::fs::create_dir(&dante_output_dir).expect("Cannot create directory.");
    }

    let task_annotation = Task::perform(
        async_tasks::run_annotation(motif_file.to_path_buf(), bam_file.to_path_buf(), output_file.clone()),
        Message::AnalysisProgress
    );
    let task_genotyping = Task::perform(
        async_tasks::run_genotyping(output_file, dante_output_dir),
        Message::AnalysisProgress
    );

    data.message_line = "Analysis started. It might take some time.".to_string();

    return task_annotation.chain(task_genotyping);
}

fn open_results(state: &mut Data) {
    let output = state.analysis_path.to_string_lossy().to_string() + "/report.html";
    opener::open(output).unwrap();
}

fn draw_open_button<'a>(state: &Data) -> Element<'a, Message> {
    let output = state.analysis_path.to_string_lossy().to_string() + "/report.html";
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

fn view_report<'a>(mut content: Column<'a, Message>, data: &'a Data, size: Size) -> Column<'a, Message> {
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

    const PAD2: Padding = Padding { bottom: 0.0, top: 10.0, right: 25.0, left: 0.0 };
    const PAD3: Padding = Padding { bottom: 0.0, top: 10.0, right: 5.0, left: 0.0 };
    let report_types = [ReportType::OnePage, ReportType::Summary, ReportType::Result, ReportType::Technical];

    let r = row![
        container(text("")).width(160),
        container(button("View").on_press(Message::EditResults(data.clone()))).padding(PAD2),
        container(button("Print").on_press(Message::Print)).padding(PAD3),
        container(pick_list(report_types, data.selected_report, Message::SetReport).placeholder("report type")).padding(PAD3),
        container(text(data.message_line2.clone())).padding(PAD2),
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

fn print_report(data: &mut Data) {
    let Ok(Some(output_pdf)) = FileDialog::new().show_save_single_file() else { return; };
    data.message_line2 = format!(
        "{} report saved to {}", data.selected_report.unwrap(), output_pdf.to_string_lossy()
    );

    // println!("{:#?}", data);
    // println!("{:?}", output_pdf);
    let typst_template = include_str!("../assets/templates/report_onepage.typ");

    use typst_as_lib::package_resolver::FileSystemCache;
    use typst_as_lib::package_resolver::PackageResolver;
    use typst_as_lib::typst_kit_options::TypstKitFontOptions;
    use typst_as_lib::TypstEngine;

    let typst_cache = App::DATA_DIR.to_string() + "/typst_cache";
    let pkg_resolver = PackageResolver::builder().cache(FileSystemCache(PathBuf::from(typst_cache))).build();

    let template = TypstEngine::builder()
        .main_file(typst_template)
        .search_fonts_with(TypstKitFontOptions::default())
        .add_file_resolver(pkg_resolver)
        .with_file_system_resolver(App::DATA_DIR)
        .build();

    let doc = template.compile().output.expect("typst::compile() returned an error!");

    let options = Default::default();

    let pdf = typst_pdf::pdf(&doc, &options).expect("Could not generate pdf.");
    std::fs::write(output_pdf, pdf).expect("Could not write pdf.");
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(super) enum ReportType {
    OnePage,
    Summary,
    Result,
    Technical,
}

impl std::fmt::Display for ReportType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::OnePage => "One-page",
            Self::Summary => "Summary",
            Self::Result => "Result",
            Self::Technical => "Technical",
        })
    }
}
