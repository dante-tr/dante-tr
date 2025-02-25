use std::path::PathBuf;
use native_dialog::FileDialog;

use iced::alignment::{Horizontal, Vertical};
use iced::widget::{
    button, checkbox, column, container, horizontal_rule, horizontal_space, pick_list,
    row, text, text_input, scrollable,
    // vertical_space, tooltip,
};
use iced::widget::{Column, Row};
use iced::{Element, Length};

use crate::{App, ContentPage};

#[derive(Debug, Clone)]
pub(super) enum Message {
    Back,
    SetMotifs(MotifFile),

    ProbandSetSex(Sex),
    ProbandSelect,
    ProbandEdit(String),

    RelativeAdd,
    RelativeRemove(usize),
    RelativeSetRelation(usize, Relation),
    RelativeSelect(usize),
    RelativeEdit(usize, String),

    Analyze,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct Data {
    pub path: PathBuf,
    pub analysis_name: String,
    pub selected: Option<MotifFile>,
    pub selected_file: Option<PathBuf>,
    pub proband_bam: Option<PathBuf>,
    pub proband_sex: Option<Sex>,

    pub relatives: Vec<(Option<PathBuf>, Option<Relation>)>,
}

impl Data {
    pub(super) fn init(path: PathBuf, analysis_name: String) -> ContentPage {
        ContentPage::AnalysisFamily(Data {
            path, analysis_name,
            selected: None, selected_file: None,
            proband_bam: None, proband_sex: None,
            relatives: vec![(None, None)]
        })
    }

    pub(super) fn view(&self) -> Element<Message> {
        println!("{:?}", self);
        let mut content = column![].align_x(Horizontal::Center).width(Length::Fill).padding(App::PAD1);

        content = make_header(content);
        content = make_form(content, self);
        content = content.push(horizontal_rule(2));
        content = make_report(content, self);

        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return scrollable(content).into();
    }

    pub(super) fn update(&mut self, m: Message) {
        // println!("{:?}\n{:?}", self, m);
        match m {
            Message::SetMotifs(motif_file) => { update_motif_selection(self, motif_file); },
            Message::ProbandSetSex(role) => { self.proband_sex = Some(role) },
            Message::ProbandSelect => { 
                if let Ok(Some(x)) = FileDialog::new().show_open_single_file() {
                    self.proband_bam = Some(x);
                }
            }
            Message::ProbandEdit(text) => { println!("{}", text); }

            Message::RelativeSetRelation(idx, role) => { self.relatives[idx].1 = Some(role); },
            Message::RelativeAdd => { self.relatives.push((None, None)); }
            Message::RelativeRemove(idx) => { self.relatives.remove(idx); }
            Message::RelativeSelect(idx) => {
                if let Ok(Some(x)) = FileDialog::new().show_open_single_file() {
                    self.relatives[idx].0 = Some(x);
                }
            }
            Message::RelativeEdit(idx, text) => { println!("{} {}", idx, text); }

            Message::Analyze => { println!("Analyze"); }

            Message::Back => { unreachable!() /* implemented in App::update */ }
        }
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
    // content = content.push(vertical_space());
    return content;
}

fn make_analyze_button(data: &Data) -> Element<Message> {
    let inactive = row![
        container(text("")).width(160), container(button("Analyze")), horizontal_space(),
    ].padding(10).align_y(Vertical::Center).into();

    let active = row![
        container(text("")).width(160),
        container(button("Analyze").on_press(Message::Analyze)),
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
        Relation::Mother, Relation::Father,
        Relation::Sister, Relation::Brother,
        Relation::Daughter, Relation::Son,
        Relation::Mate
    ];

    let (path, relation) = &data.relatives[0];
    let filename = path.clone().unwrap_or_default().to_string_lossy().to_string();
    let text_message = |x| { Message::RelativeEdit(0, x) };
    let pick_message = |x| { Message::RelativeSetRelation(0, x) };
    let first_row = row![
        container(text("Relatives: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text_input("Type path or click search...", &filename).on_input(text_message)).padding(App::PAD1),
        container(pick_list(choices, *relation, pick_message).placeholder("relation")).padding(App::PAD1),
        container(button("Search").on_press(Message::RelativeSelect(0))).padding(App::PAD1)
    ].padding(10).align_y(Vertical::Center);

    content = content.push(first_row);

    for (i, (path, relation)) in data.relatives.iter().enumerate().skip(1) {
        let filename = path.clone().unwrap_or_default().to_string_lossy().to_string();
        let text_message = move |x| { Message::RelativeEdit(i, x) };
        let pick_message = move |x| { Message::RelativeSetRelation(i, x) };
        let next_row = row![
            container(button("Remove").on_press(Message::RelativeRemove(i))).width(160).align_x(Horizontal::Right).padding(App::PAD1),
            container(text_input("Type path or click search...", &filename).on_input(text_message)).padding(App::PAD1),
            container(pick_list(choices, *relation, pick_message).placeholder("relation")).padding(App::PAD1),
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
    let sex = [ Sex::Male, Sex::Female, Sex::Intersex ];

    row![
        container(text("Proband: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text_input("Type path or click search...", &proband).on_input(Message::ProbandEdit)).padding(App::PAD1),
        container(pick_list(sex, data.proband_sex, Message::ProbandSetSex).placeholder("sex")).padding(App::PAD1),
        container(button("Search").on_press(Message::ProbandSelect)).padding(App::PAD1)
    ].padding(10).align_y(Vertical::Center)
}

fn make_motif_selection(selected: Option<MotifFile>, selected_file: &Option<PathBuf>) -> Row<Message> {
    let motif_files = [
        MotifFile::STRSet_20220902,
        MotifFile::Custom
    ];

    let content = match selected {
        None => {
            row![
                container(text("Motifs: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
                container(pick_list(motif_files, selected, Message::SetMotifs).placeholder("type")),
                container(text("")).width(Length::Fill).align_x(Horizontal::Left)
            ].padding(10).align_y(Vertical::Center)
        }
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
            /* TODO: store some path */
            // App::DATA_DIR + "/somepath/motifs.tsv"
            data.selected_file = Some(PathBuf::from(""));
            data.selected = Some(motif_file);
        }
        MotifFile::Custom => {
            if let Ok(Some(x)) = FileDialog::new().show_open_single_file() {
                data.selected_file = Some(x);
                data.selected = Some(motif_file);
            }
        }
    }
}

fn make_report<'a>(mut content: Column<'a, Message>, data: &'a Data) -> Column<'a, Message> {
    // content = content.push(vertical_space());
    let tmp = false;
    content = content.push(row![
        container(text("Filter: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        checkbox("All", tmp),
        checkbox("ALS", tmp), checkbox("DM2", tmp), checkbox("OPDM", tmp),
        checkbox("ALS", tmp), checkbox("ALS", tmp), checkbox("ALS", tmp),
        checkbox("ALS", tmp), checkbox("ALS", tmp), checkbox("ALS", tmp),
        checkbox("ALS", tmp), checkbox("ALS", tmp), checkbox("ALS", tmp),
        checkbox("ALS", tmp),
        horizontal_space(),
    ].padding(10).align_y(Vertical::Center));

    content = content.push(row![
        container(text("")).width(160).align_x(Horizontal::Right),
        checkbox("All", tmp),
        checkbox("ALS", tmp), checkbox("DM2", tmp), checkbox("OPDM", tmp),
        checkbox("ALS", tmp), checkbox("ALS", tmp), checkbox("ALS", tmp),
        checkbox("ALS", tmp), checkbox("ALS", tmp), checkbox("ALS", tmp),
        checkbox("ALS", tmp), checkbox("ALS", tmp), checkbox("ALS", tmp),
        checkbox("ALS", tmp),
        horizontal_space(),
    ].padding(10).align_y(Vertical::Center));

    content = content.push(row![
        container(text("")).width(160),
        container(button("View")),
        container(button("Print")).padding(10),
        horizontal_space(),
    ].padding(10).align_y(Vertical::Center));
    // content = content.push(vertical_space());
    return content;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum MotifFile {
    #[allow(non_camel_case_types)]
    STRSet_20220902,
    Custom,
}

impl std::fmt::Display for MotifFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::STRSet_20220902 => "STRSet_20220902",
            Self::Custom => "custom",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Sex {
    Male,
    Female,
    Intersex,
}

impl std::fmt::Display for Sex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Male => "male",
            Self::Female => "female",
            Self::Intersex => "intersex"
        })
    }
}


