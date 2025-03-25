use native_dialog::FileDialog;
use serde::{Serialize, Deserialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, checkbox, pick_list, text, text_input, tooltip};
use iced::widget::{column, container, horizontal_rule, horizontal_space, row, scrollable};
use iced::widget::{Column, Row};
use iced::{Element, Length, Padding, Size, Task};

use crate::async_tasks;
use crate::pdf_reporting;
use crate::{App, ContentPage, MotifFile, Sex};

const PSIZE: u16 = 100;

#[derive(Debug, Clone)]
pub(super) enum Message {
    Back,
    SetMotifs(MotifFile),
    MotifCheckbox(usize, bool),
    MotifGroupbox(usize, bool),

    ProbandSetSex(Sex),
    ProbandSelect,
    ProbandEdit(String),

    RelativeAdd,
    RelativeRemove(usize),
    RelativeSetRelation(usize, Relation),
    RelativeSelect(usize),
    RelativeEdit(usize, String),

    Analyze,
    AnalysisProgress(String),
    Print,
}

#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub(super) struct Data {
    path: PathBuf,
    analysis_name: String,
    selected: Option<MotifFile>,
    selected_file: Option<PathBuf>,

    proband_bam: Option<PathBuf>,
    proband_sex: Option<Sex>,

    relatives: Vec<(Option<PathBuf>, Option<Relation>)>,

    motifs: Vec<(bool, String, Vec<String>, String)>,
    groups: Vec<(bool, String)>,

    progress: String,
    tasks_done: usize,
    tasks_todo: usize,
}

impl Data {
    pub(super) fn view(&self, size: Size) -> Element<Message> {
        let mut content = column![].align_x(Horizontal::Center).width(Length::Fill).padding(App::PAD1);

        content = make_header(content);
        content = make_form(content, self);
        content = content.push(container(horizontal_rule(2)).padding(25));
        content = make_report(content, self, size);

        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return scrollable(content).into();
    }

    #[rustfmt::skip]
    pub(super) fn update(&mut self, m: Message) -> Task<Message> {
        match m {
            Message::SetMotifs(motif_file)
                => { update_motif_selection(self, motif_file); Task::none() },
            Message::ProbandSetSex(role)
                => { self.proband_sex = Some(role); Task::none() },
            Message::ProbandSelect
                => { select_file(&mut self.proband_bam); Task::none() }
            Message::ProbandEdit(text)
                => { println!("{}", text); todo!() }
            Message::RelativeSetRelation(idx, role)
                => { self.relatives[idx].1 = Some(role); Task::none() },
            Message::RelativeAdd
                => { self.relatives.push((None, None)); Task::none() }
            Message::RelativeRemove(idx)
                => { self.relatives.remove(idx); Task::none() }
            Message::RelativeSelect(idx)
                => { select_file(&mut self.relatives[idx].0); Task::none() }
            Message::RelativeEdit(idx, text)
                => { println!("{} {}", idx, text); todo!() }
            Message::MotifCheckbox(idx, checked)
                => { self.motifs[idx].0 = checked; Task::none() }
            Message::MotifGroupbox(idx, checked)
                => { toggle_group(self, idx, checked); Task::none() }
            Message::Print
                => { print_report(self); Task::none() }
            Message::Back
                => { unreachable!() /* implemented in App::update */ }

            Message::Analyze => { analyze(self) }
            Message::AnalysisProgress(p) => {
                self.tasks_done += 1;
                if self.tasks_done == self.tasks_todo {
                    self.progress = format!(
                        "{}/{} tasks done. Analysis finished successfully.",
                        self.tasks_done, self.tasks_todo
                    )
                } else {
                    self.progress = format!(
                        "{}/{} tasks done. {}",
                        self.tasks_done, self.tasks_todo, p
                    )
                };
                Task::none()
            }
        }
    }

    pub(super) fn init(path: PathBuf, analysis_name: String) -> ContentPage {
        let data = Data {
            path, analysis_name, relatives: vec![(None, None)], ..Default::default()
        };
        data.save();
        ContentPage::AnalysisFamily(data)
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

fn select_file(place: &mut Option<PathBuf>) {
    if let Ok(Some(x)) = FileDialog::new().show_open_single_file() {
        *place = Some(x);
    }
}

fn make_header(mut content: Column<Message>) -> Column<Message> {
    content = content.push(row![
        container(button("Back").on_press(Message::Back)).width(100),
        container(text("Family analysis").size(App::H1_SIZE)).align_x(Horizontal::Center).width(Length::Fill),
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

    content = make_relatives(content, data);
    content = content.push(make_analyze_button(data));
    return content;
}

fn make_analyze_button(data: &Data) -> Element<Message> {
    let inactive = row![
        container(text("")).width(160), container(button("Analyze")), horizontal_space(),
    ].padding(10).align_y(Vertical::Center).into();

    let active = row![
        container(text("")).width(160),
        container(button("Analyze").on_press(Message::Analyze)),
        container(text(data.progress.clone())).padding(App::PAD2),
        horizontal_space(),
    ].padding(10).align_y(Vertical::Center).into();

    if data.selected.is_none() { return inactive; }
    if data.proband_bam.is_none() { return inactive; }
    if data.proband_sex.is_none() { return inactive; }

    for (relation, bam) in &data.relatives {
        if bam.is_none() { return inactive; }
        if relation.is_none() { return inactive; }
    }
    return active;
}

fn make_relatives<'a>(mut content: Column<'a, Message>, data: &'a Data) -> Column<'a, Message> {
    let choices = [
        Relation::Mother, Relation::Father, Relation::Sister, Relation::Brother,
        Relation::Daughter, Relation::Son, Relation::Mate
    ];

    let (path, relation) = &data.relatives[0];
    let filename = path.clone().unwrap_or_default().to_string_lossy().to_string();
    let text_message = |x| Message::RelativeEdit(0, x);
    let pick_message = |x| Message::RelativeSetRelation(0, x);
    let first_row = row![
        container(text("Relatives: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text_input("Type path or click search...", &filename).on_input(text_message)).padding(App::PAD1),
        container(pick_list(choices, *relation, pick_message).placeholder("relation").width(PSIZE)).padding(App::PAD1),
        container(button("Search").on_press(Message::RelativeSelect(0))).padding(App::PAD1)
    ].padding(10).align_y(Vertical::Center);

    content = content.push(first_row);

    for (i, (path, relation)) in data.relatives.iter().enumerate().skip(1) {
        let filename = path.clone().unwrap_or_default().to_string_lossy().to_string();
        let text_message = move |x| Message::RelativeEdit(i, x);
        let pick_message = move |x| Message::RelativeSetRelation(i, x);
        let next_row = row![
            container(button("Remove").on_press(Message::RelativeRemove(i))).width(160).align_x(Horizontal::Right).padding(App::PAD1),
            container(text_input("Type path or click search...", &filename).on_input(text_message)).padding(App::PAD1),
            container(pick_list(choices, *relation, pick_message).placeholder("relation").width(PSIZE)).padding(App::PAD1),
            container(button("Search").on_press(Message::RelativeSelect(i))).padding(App::PAD1)
        ].padding(10).align_y(Vertical::Center);

        content = content.push(next_row);
    }

    content = content.push(row![
        container(text("")).width(160).align_x(Horizontal::Right),
        container(button("Add relative").on_press(Message::RelativeAdd)).padding(App::PAD1),
        container(text("")).width(Length::Fill).align_x(Horizontal::Left)
    ].padding(10).align_y(Vertical::Center));
    return content;
}

fn make_proband_row(data: &Data) -> Row<Message> {
    let proband = data.proband_bam.clone().unwrap_or_default().to_string_lossy().to_string();
    let sex = [Sex::Male, Sex::Female, Sex::Unknown];

    row![
        container(text("Proband: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text_input("Type path or click search...", &proband).on_input(Message::ProbandEdit)).padding(App::PAD1),
        container(pick_list(sex, data.proband_sex, Message::ProbandSetSex).placeholder("sex").width(PSIZE)).padding(App::PAD1),
        container(button("Search").on_press(Message::ProbandSelect)).padding(App::PAD1)
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

fn update_motif_selection(data: &mut Data, motif_file: MotifFile) {
    match motif_file {
        MotifFile::STRSet_20220902 => {
            let path = PathBuf::from(App::DATA_DIR.to_string() + "/STRSet_20220902.tsv");
            data.selected = Some(motif_file);
            data.selected_file = Some(path);
            data.motifs = parse_motifs(data.selected_file.as_ref().unwrap());
            data.groups = get_groups(data.motifs.as_ref());
        },
        MotifFile::STRSet_20250311 => {
            let path = PathBuf::from(App::DATA_DIR.to_string() + "/STRSet_20250311.tsv");
            data.selected = Some(motif_file);
            data.selected_file = Some(path);
            data.motifs = parse_motifs(data.selected_file.as_ref().unwrap());
            data.groups = get_groups(data.motifs.as_ref());
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

fn toggle_group(data: &mut Data, idx: usize, checked: bool) {
    data.groups[idx].0 = checked;
    let group = data.groups[idx].1.clone();
    for x in &mut data.motifs { if x.2.contains(&group) { x.0 = checked; } }
}

fn parse_motifs(path: &Path) -> Vec<(bool, String, Vec<String>, String)> {
    let file = File::open(path).expect("Cannot find motif file.");
    let reader = BufReader::new(file);

    let mut result = Vec::new();
    for line in reader.lines() {
        let line = line.expect("Cannot read line from motif file.").trim().to_string();
        let split: Vec<_> = line.split('\t').collect();

        let id = split[0].to_string();
        // let hgvs = split[1];
        let groups = split[2].split(',').map(|x| x.to_string()).collect();
        let description = split[3].to_string();
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

fn analyze(data: &mut Data) -> Task<Message> {
    data.save();

    let mut tasks = Vec::new();
    let Some(ref motif_file) = data.selected_file else { unreachable!() }; 
    let Some(ref bam_file) = data.proband_bam else { unreachable!() };

    let task_proband = get_chain(motif_file, bam_file, data.path.clone());
    tasks.push(task_proband);

    for relative in &data.relatives {
        let Some(ref bam_file) = relative.0 else {unreachable!() };
        let task_chain = get_chain(motif_file, bam_file, data.path.clone());
        tasks.push(task_chain);

    }
    data.tasks_done = 0;
    data.tasks_todo = tasks.len() * 2;
    data.progress = format!(
        "{}/{} tasks done. Analysis started. It might take some time.",
        data.tasks_done, data.tasks_todo
    );

    return Task::batch(tasks);
}

fn get_chain(motif_file: &Path, bam_file: &Path, mut output_file: PathBuf) -> Task<Message> {
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
    let task_chain = task_annotation.chain(task_genotyping);
    return task_chain;
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(super) enum Relation {
    Mother,
    Father,
    Sister,
    Brother,
    Daughter,
    Son,
    Mate,
}

impl std::fmt::Display for Relation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Mother => "mother",
            Self::Father => "father",
            Self::Sister => "sister",
            Self::Brother => "brother",
            Self::Daughter => "daughter",
            Self::Son => "son",
            Self::Mate => "mate",
        })
    }
}

fn print_report(data: &mut Data) {
    pdf_reporting::simple_report(data);
    opener::open("typst_report.pdf").unwrap();
    println!("{:?}", data);
    
}
