use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, column, container, horizontal_rule, row, text, text_input};
use iced::widget::pick_list;
use iced::Element;

use crate::App;
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
    println!("{:?}", state);
    println!("{:?}", analysis);

    use ContentPage as CP;
    if let CP::WelcomeScreen(data) = &mut state.content_page {
        data.selected = Some(analysis);
    }
}

pub fn analysis_create(state: &mut App) {
    use ContentPage as CP;
    let CP::WelcomeScreen(data) = state.content_page else { unreachable!() };
    let Some(x) = data.selected else { unreachable!() };
    match x {
        Analysis::Single => { state.content_page = CP::AnalysisSingle; },
        Analysis::Family => { state.content_page = CP::AnalysisFamily; }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct Data {
    selected: Option<Analysis>,
}

pub fn view(_state: &App, data: Data) -> Element<Message> {
    use App as S;
    let analyses = [
        Analysis::Single,
        Analysis::Family,
    ];

    let dropdown1 = pick_list(analyses, data.selected, Message::AnalysisSelected).placeholder("type");
    let button1 = match data.selected {
        None => { button("Create") },
        Some(_) => { button("Create").on_press(Message::AnalysisCreate) }
    };

    let value = String::new();
    // text_input("Type path or click search...", &filename_str).on_input(on_input).font(App::BOLD_MONO),

    column![
        row![
            container(text("Create new analysis: ").width(160).align_x(Horizontal::Right)).padding(S::PAD1),
            container(text_input("name", &value).font(App::BOLD_MONO)).padding(S::PAD1),
            container(dropdown1).padding(S::PAD1),
            container(button1).padding(S::PAD2),
        ].padding(10.0).align_y(Vertical::Center),
        horizontal_rule(2),
        row![
            container(text("Previous analyses: ").width(160).align_x(Horizontal::Right)).padding(App::PAD1)
        ].padding(10.0).align_y(Vertical::Center),

    ].width(720.0).align_x(Horizontal::Left).into()
}


