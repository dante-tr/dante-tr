use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, column, container, horizontal_rule, row, text};
use iced::widget::pick_list;
use iced::Element;

use crate::State;
use crate::Message;
use crate::NewState;

pub fn view(state: &State) -> Element<Message> {
    use State as S;
    let analyses = [
        NewState::AnalysisSingle,
        NewState::AnalysisFamily,
    ];

    let dropdown1 = pick_list(analyses, state.analysis, Message::AnalysisSelected).placeholder("type");

    column![
        row![
            container(text("Create new analysis: ").width(180).align_x(Horizontal::Right)).padding(S::PAD1),
            container(dropdown1).padding(S::PAD1),
            container(button("Create")).padding(S::PAD2),
            //.on_press(on_press)).padding(State::PAD2),
        ].padding(10.0).align_y(Vertical::Center),
        horizontal_rule(2),
        row![
            container(text("Previous analyses:").width(180).align_x(Horizontal::Right)).padding(State::PAD1)
        ].padding(10.0).align_y(Vertical::Center),

    ].width(720.0).align_x(Horizontal::Left).into()
}


