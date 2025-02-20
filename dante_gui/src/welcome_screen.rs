use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, column, container, horizontal_rule, row, text, text_input, tooltip};
use iced::widget::pick_list;
use iced::Element;
use std::fs;
use std::path::{Path, PathBuf};
use std::ffi::OsStr;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::{analysis_single, App, ContentPage};
use crate::analysis_family::Data as FamilyData;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Data {
    name: String,
    selected: Option<Analysis>,
}


#[derive(Debug, Clone)]
pub(crate) enum Message {
    AnalysisNamed(String),
    AnalysisSelected(Analysis),
    CreateAnalysis,
    AnalysisReopen(PathBuf),
    AnalysisDelete(PathBuf),
}

pub(crate) fn view(data: & Data) -> Element<Message> {
    let analyses = [
        Analysis::Single,
        Analysis::Family,
    ];

    let dropdown1 = pick_list(analyses, data.selected, Message::AnalysisSelected).placeholder("type");
    let button1 = make_button1(data);
    let previous_analyses = make_previous_list(data);

    column![
        row![
            container(text("Create new analysis: ").width(160).align_x(Horizontal::Right)).padding(App::PAD1),
            container(text_input("name", &data.name).on_input(Message::AnalysisNamed)).padding(App::PAD1),
            container(dropdown1).padding(App::PAD1),
            container(button1).padding(App::PAD2),
        ].padding(10.0).align_y(Vertical::Center),
        horizontal_rule(2),
        previous_analyses
    ].align_x(Horizontal::Center).into()
}

pub(crate) fn update(data: &mut Data, m: Message) {
    match m {
        Message::AnalysisSelected(analysis) => { data.selected = Some(analysis); },
        Message::AnalysisNamed(name) => { data.name = name; },
        Message::AnalysisDelete(path) => { fs::remove_dir_all(path).unwrap(); }
        Message::CreateAnalysis => { unreachable!() /* implemented in App::update */ },
        Message::AnalysisReopen(_) => { unreachable!() /* implemented in App::update */ },
    }
}

pub(crate) fn analysis_create(state: &mut App) {
    use ContentPage as CP;
    let (atype, name) = match &state.content_page {
        CP::WelcomeScreen(Data{selected: Some(atype), name}) => { (atype, name) },
        _ => { unreachable!() }
    };
    let name = name.to_string();

    // TODO: make it human readable? Or not?
    let time = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Clock may have gone backwards. What are you doing?")
        .as_secs().to_string();
    let path: PathBuf = format!("{}/analyses/{}_{}_{}", App::DATA_DIR, time, name, atype).into();
    mkdir_p(&path);

    match atype {
        Analysis::Single => { 
            state.content_page = CP::AnalysisSingle(analysis_single::Data {
                analysis_name: name,
                ..Default::default()
            });
        },
        Analysis::Family => { state.content_page = FamilyData::init(path, name); }
    }
}

pub(crate) fn analysis_reopen(state: &mut App, path: PathBuf) {
    let (_, name, atype) = parse_analysis_dir(&path);
    use ContentPage as CP;
    match atype.as_str() {
        "single" => {
            state.content_page = CP::AnalysisSingle(analysis_single::Data {
                analysis_name: name,
                ..Default::default()
            })
        },
        "family" => { state.content_page = FamilyData::init(path, name); },
        _ => { unreachable!() }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Analysis {
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
    const TYPE_WIDTH: u16 = 150;
    const ACTION_WIDTH: u16 = 150;
    let mut result = column![].align_x(Horizontal::Center).width(NAME_WIDTH + TYPE_WIDTH + ACTION_WIDTH + 10);
    result = result.push(
        row![
            container(text("Previous analyses").size(22)).align_x(Horizontal::Center)
        ].padding(15.0).align_y(Vertical::Center)
    );

    result = result.push(
        row![
            container(text("Analysis name").width(NAME_WIDTH).align_x(Horizontal::Center)).padding(App::PAD1),
            container(text("type").width(TYPE_WIDTH).align_x(Horizontal::Center)).padding(App::PAD1),
            container(text("actions").width(ACTION_WIDTH).align_x(Horizontal::Center)),
        ].padding(5.0).align_y(Vertical::Center)
    );
    result = result.push(horizontal_rule(2));

    for path in paths.into_iter().rev() {
        let (_, analysis_name, analysis_type) = parse_analysis_dir(&path);
        let path_name = path.clone().into_os_string().into_string().unwrap();

        result = result.push(
            row![
                tooltip(
                    container(text(analysis_name).width(NAME_WIDTH).align_x(Horizontal::Center)).padding(App::PAD1),
                    text(path_name),
                    tooltip::Position::FollowCursor
                ),
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
    if data.selected.is_none() { return button("Create").into() }
    if data.name.is_empty() { return button("Create").into() }
    return button("Create").on_press(Message::CreateAnalysis).into()
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
