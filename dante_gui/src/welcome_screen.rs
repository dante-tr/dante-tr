use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, column, container, horizontal_rule, row, text, text_input, tooltip};
use iced::widget::pick_list;
use iced::Element;
use std::fs;
use std::path::PathBuf;

use crate::{mkdir_p, App};
use crate::Message;
use crate::ContentPage;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Analysis {
    Single,
    Family,
}

impl std::fmt::Display for Analysis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Single => "Single",
            Self::Family => "Family",
        })
    }
}

pub fn analysis_set(state: &mut App, analysis: Analysis) {
    use ContentPage as CP;
    if let CP::WelcomeScreen(data) = &mut state.content_page {
        data.selected = Some(analysis);
    }
}

pub fn analysis_create(state: &mut App) {
    use ContentPage as CP;
    let CP::WelcomeScreen(ref data) = state.content_page else { unreachable!() };
    let Some(x) = data.selected else { unreachable!() };

    let path = App::DATA_DIR.to_string() + "/analyses/" + &data.name;
    mkdir_p(path);

    match x {
        Analysis::Single => { state.content_page = CP::AnalysisSingle; },
        Analysis::Family => { state.content_page = CP::AnalysisFamily; }
    }
}

pub fn analysis_name(state: &mut App, name: String) {
    use ContentPage as CP;
    let CP::WelcomeScreen(ref mut data) = state.content_page else { unreachable!() };
    data.name = name
}

pub(crate) fn analysis_reopen(_state: &mut App, path: PathBuf) {

    println!("{:?}", path);
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Data {
    name: String,
    selected: Option<Analysis>,
}

pub fn view<'a>(state: &'a App, data: &'a Data) -> Element<'a, Message> {
    use App as S;
    let analyses = [
        Analysis::Single,
        Analysis::Family,
    ];

    let dropdown1 = pick_list(analyses, data.selected, Message::AnalysisSelected).placeholder("type");
    let button1 = make_button1(data);
    let previous_analyses = previous(state);

    column![
        row![
            container(text("Create new analysis: ").width(160).align_x(Horizontal::Right)).padding(S::PAD1),
            container(text_input("name", &data.name).on_input(Message::AnalysisNamed).font(App::BOLD_MONO)).padding(S::PAD1),
            container(dropdown1).padding(S::PAD1),
            container(button1).padding(S::PAD2),
        ].padding(10.0).align_y(Vertical::Center),
        horizontal_rule(2),
        previous_analyses
    ].width(720.0).align_x(Horizontal::Left).into()
}

fn make_button1(data: &Data) -> Element<Message> {
    if data.selected.is_none() { return button("Create").into() }
    if data.name.is_empty() { return button("Create").into() }
    return button("Create").on_press(Message::AnalysisCreate).into()
}

fn previous(_state: &App) -> Element<Message> {
    let mut result = column![];

    let analyses_dir = PathBuf::from(App::DATA_DIR.to_string() + "/analyses/");
    let paths: Vec<PathBuf> = if analyses_dir.exists() {
        fs::read_dir(analyses_dir).unwrap().map(|x| x.unwrap().path()).collect()
    } else {
        Vec::new()
    };

    result = result.push(
        row![
            container(text("Previous analyses: ").width(160).align_x(Horizontal::Right)).padding(App::PAD1)
        ].padding(10.0).align_y(Vertical::Center)
    );
    for path in paths {
        let a: String = path.file_name().unwrap().to_owned().into_string().unwrap();
        let b: String = String::from("Single");
        let c: String = path.clone().into_os_string().into_string().unwrap();
        let y = row![
            container(button("Load").on_press(Message::AnalysisReopen(path))).width(150).align_x(Horizontal::Center).padding(App::PAD1),
            tooltip(
                text(a).width(100).align_x(Horizontal::Left),
                container(text(c)).padding(App::PAD1).style(container::rounded_box),
                tooltip::Position::FollowCursor
            ),
            container(text(b).width(100).align_x(Horizontal::Left)).padding(App::PAD1),
            // container(text(c).width(500).align_x(Horizontal::Left)).padding(App::PAD1),
        ].padding(5.0).align_y(Vertical::Center);
        result = result.push(y);
    }
    return result.into();
}

