use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, column, container, horizontal_rule, row, text, text_input, tooltip};
use iced::widget::pick_list;
use iced::Element;
use std::fs;
use std::path::{Path, PathBuf};
use std::ffi::OsStr;
use std::time::SystemTime;
use chrono::{DateTime, Local};

use crate::{App, ContentPage};
use crate::analysis_family::Data as FamilyData;
use crate::analysis_single::Data as SingleData;

#[derive(Debug, Clone)]
pub(super) enum Message {
    AnalysisNamed(String),
    AnalysisSelected(Analysis),
    CreateAnalysis(String, Analysis),
    AnalysisReopen(PathBuf),
    AnalysisDelete(PathBuf),
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct Data {
    name: String,
    selected: Option<Analysis>,
}

impl Data {
    pub(super) fn view(&self) -> Element<Message> {
        let mut content = column![].align_x(Horizontal::Center);

        let analyses = [
            Analysis::Single,
            Analysis::Family,
        ];

        let dropdown1 = pick_list(analyses, self.selected, Message::AnalysisSelected).placeholder("type");
        let button1 = make_button1(self);
        let previous_analyses = make_previous_list(self);

        content = content.push(
            row![
                container("").width(100),
                container(text("Create new analysis: ").width(160).align_x(Horizontal::Right)).padding(App::PAD1),
                container(text_input("name", &self.name).on_input(Message::AnalysisNamed)).padding(App::PAD1),
                container(dropdown1).padding(App::PAD1),
                container(button1).padding(App::PAD2),
                container("").width(100),
            ].padding(25).align_y(Vertical::Center)
        );
        content = content.push(horizontal_rule(2));
        content = content.push(previous_analyses);

        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return content.into();
    }

    pub(super) fn update(&mut self, m: Message) {
        match m {
            Message::AnalysisSelected(analysis) => { self.selected = Some(analysis); },
            Message::AnalysisNamed(name) => { self.name = name; },
            Message::AnalysisDelete(path) => { fs::remove_dir_all(path).unwrap(); }
            Message::CreateAnalysis(_, _) => { unreachable!() /* implemented in App::update */ }
            Message::AnalysisReopen(_) => { unreachable!() /* implemented in App::update */ },
        }
    }
}

pub(super) fn analysis_create(name: String, atype: Analysis) -> ContentPage {
    let time: DateTime<Local> = SystemTime::now().into();
    let time = time.format("%Y-%m-%d-%H-%M-%S");

    let path: PathBuf = format!("{}/analyses/{}_{}_{}", App::DATA_DIR, time, name, atype).into();
    mkdir_p(&path);

    match atype {
        Analysis::Single => { return SingleData::init(path, name); },
        Analysis::Family => { return FamilyData::init(path, name); }
    }
}

pub(super) fn analysis_reopen(path: PathBuf) -> ContentPage {
    let (_, _, atype) = parse_analysis_dir(&path);
    use ContentPage as CP;
    match atype.as_str() {
        "single" => { return CP::AnalysisSingle(SingleData::load(path)); },
        "family" => { return CP::AnalysisFamily(FamilyData::load(path)); },
        _ => { unreachable!() }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Analysis {
    Single,
    Family,
}

impl std::fmt::Display for Analysis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Single => "single",
            Self::Family => "family",
        })
    }
}

fn make_previous_list(_: &Data) -> Element<Message> {
    let analyses_dir = PathBuf::from(App::DATA_DIR.to_string() + "/analyses/");
    let mut paths: Vec<PathBuf> = if analyses_dir.exists() {
        fs::read_dir(analyses_dir).unwrap().map(|x| x.unwrap().path()).collect()
    } else {
        Vec::new()
    };
    paths.sort();

    const NAME_WIDTH: u16 = 350;
    const TIME_WIDTH: u16 = 200;
    const TYPE_WIDTH: u16 = 150;
    const ACTION_WIDTH: u16 = 150;
    let mut result = column![].align_x(Horizontal::Center).width(NAME_WIDTH + TIME_WIDTH + TYPE_WIDTH + ACTION_WIDTH + 10);
    result = result.push(
        row![
            container(text("Previous analyses").size(22)).align_x(Horizontal::Center)
        ].padding(15.0).align_y(Vertical::Center)
    );

    result = result.push(
        row![
            container(text("Analysis name").width(NAME_WIDTH).align_x(Horizontal::Center)).padding(App::PAD1),
            container(text("time").width(TIME_WIDTH).align_x(Horizontal::Center)).padding(App::PAD1),
            container(text("type").width(TYPE_WIDTH).align_x(Horizontal::Center)).padding(App::PAD1),
            container(text("actions").width(ACTION_WIDTH).align_x(Horizontal::Center)),
        ].padding(5.0).align_y(Vertical::Center)
    );
    result = result.push(horizontal_rule(2));

    for path in paths.into_iter().rev() {
        let (analysis_time, analysis_name, analysis_type) = parse_analysis_dir(&path);
        let path_name = path.clone().into_os_string().into_string().unwrap();

        let parts: Vec<_> = analysis_time.split("-").collect();
        let analysis_time = format!("{}-{}-{} {}:{}:{}", parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]);
        result = result.push(
            row![
                tooltip(
                    container(text(analysis_name).width(NAME_WIDTH).align_x(Horizontal::Center)).padding(App::PAD1),
                    container(text(format!("Located in: {}", path_name))).padding(5).style(container::rounded_box),
                    tooltip::Position::FollowCursor
                ),
                container(text(analysis_time).width(TIME_WIDTH).align_x(Horizontal::Center)).padding(App::PAD1),
                container(text(analysis_type).width(TYPE_WIDTH).align_x(Horizontal::Center)).padding(App::PAD1),
                column![row![
                    container(button("Load").on_press(Message::AnalysisReopen(path.clone()))).padding(App::PAD2),
                    container(button("Delete").on_press(Message::AnalysisDelete(path.clone()))).padding(App::PAD2),
                ]].width(ACTION_WIDTH).align_x(Horizontal::Center)
            ].padding(5.0).align_y(Vertical::Center)
        );
    }

    // let result = std::convert::Into::<Element<Message>>::into(result).explain(iced::Color::BLACK);
    let result = result.into();
    return result;
}

fn make_button1(data: &Data) -> Element<Message> {
    if data.name.is_empty() { return button("Create").into() }
    let msg = match (&data.name, data.selected) {
        (_, None) => { return button("Create").into(); },
        (x, Some(y)) => { Message::CreateAnalysis(x.clone(), y) }
    };
    return button("Create").on_press(msg).into()
}

fn mkdir_p<P>(dir: P)
where 
     P: AsRef<OsStr> + AsRef<Path>
{
    let path = Path::new(&dir);
    let ancestors: Vec<_> = path.ancestors().collect();
    for anc in ancestors.iter().rev().skip(1) {
        println!("{:?}", anc);
        if !anc.exists() {
            fs::create_dir(anc).expect("Cannot create directory.");
        }
    }
}

fn parse_analysis_dir(path: &Path) -> (String, String, String) {
    let filename = path.file_name().unwrap().to_owned().into_string().unwrap();
    let x: Vec<usize> = filename.match_indices("_").map(|x| x.0).collect();
    let a = x[0];
    let b = *x.last().unwrap();
    let n = filename.len();

    let analysis_time = filename[0..a].to_string();
    let analysis_name = filename[a+1..b].to_string();
    let analysis_type = filename[b+1..n].to_string();

    return (analysis_time, analysis_name, analysis_type);
}
