use std::path::PathBuf;

use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, checkbox, column, container, horizontal_rule, horizontal_space, row, text, tooltip, vertical_space};
use iced::widget::Column;
use iced::{Element, Length};

use crate::{App, ContentPage};

#[derive(Debug, Clone)]
pub(crate) enum Message {
    Back,
    Temp,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Data {
    pub path: PathBuf,
    pub analysis_name: String,
}

impl Data {
    pub(crate) fn init(path: PathBuf, analysis_name: String) -> ContentPage {
        ContentPage::AnalysisFamily(Data { path, analysis_name })
    }

    pub(crate) fn view(&self) -> Element<Message> {
        println!("{:?}", self);
        let mut content = column![].align_x(Horizontal::Center).width(Length::Fill);

        content = make_header(content);
        content = make_form(content, self);
        content = content.push(horizontal_rule(2));
        content = make_report(content, self);

        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return content.into();
    }

    pub(crate) fn update(&mut self, m: Message) {
        println!("{:?}\n{:?}", self, m)
    }
}

fn make_header(mut content: Column<Message>) -> Column<Message> {
    content = content.push(vertical_space());
    content = content.push(row![
        container(button("Back").on_press(Message::Back)).width(100),
        container(text("Family analysis").size(App::H1_SIZE)).align_x(Horizontal::Center).width(Length::Fill),
        container("").width(100),
    ].padding(10).align_y(Vertical::Center));
    content = content.push(vertical_space());
    return content;
}

fn make_form<'a>(mut content: Column<'a, Message>, data: &'a Data) -> Column<'a, Message> {
    content = content.push(vertical_space());
    content = content.push(row![
        container(text("Analysis name: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text(data.analysis_name.clone())).width(Length::Fill).align_x(Horizontal::Left)
    ].padding(10).align_y(Vertical::Center));

    content = content.push(row![
        container(text("Motifs: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text("")).width(Length::Fill).align_x(Horizontal::Left)
    ].padding(10).align_y(Vertical::Center));

    content = content.push(row![
        container(text("Proband: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text("")).width(Length::Fill).align_x(Horizontal::Left)
    ].padding(10).align_y(Vertical::Center));

    content = content.push(row![
        container(text("Relatives: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
        container(text("")).width(Length::Fill).align_x(Horizontal::Left)
    ].padding(10).align_y(Vertical::Center));

    content = content.push(row![
        container(text("")).width(160).align_x(Horizontal::Right),
        container(text("")).width(Length::Fill).align_x(Horizontal::Left)
    ].padding(10).align_y(Vertical::Center));

    content = content.push(row![
        container(text("")).width(160),
        container(button("Analyze")),
        horizontal_space(),
    ].padding(10).align_y(Vertical::Center));
    content = content.push(vertical_space());
    return content;
}

fn make_report<'a>(mut content: Column<'a, Message>, data: &'a Data) -> Column<'a, Message> {
    content = content.push(vertical_space());
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
    content = content.push(vertical_space());
    return content;
}
