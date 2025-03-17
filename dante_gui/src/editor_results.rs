use iced::alignment::{Horizontal, Vertical};
use iced::widget::{column, row, container, button, text};
use iced::widget::Column;
use iced::{Element, Size, Length};

use std::path::PathBuf;
use serde::{Serialize, Deserialize};
use serde_json::Value;

use crate::ContentPage;
use crate::App;

#[derive(Debug, Clone)]
pub(crate) enum Message {
    Exit(PathBuf),
    Save,
}

#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub(super) struct Data {
    source: PathBuf,
    motif_ids: Vec<String>,
}

impl Data {
    pub(super) fn view(&self, _size: Size) -> Element<Message> {
        let mut content = column![].align_x(Horizontal::Center);

        content = view_header(content, self.source.clone());
        content = view_results(content, self);
        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return content.into();
    }

    pub(super) fn update(&mut self, m: Message) {
        match m {
            Message::Exit(_) => { unreachable!("Implemented in App::update."); }
            Message::Save => { println!("Save"); }
        }
    }

    pub(super) fn open() -> ContentPage {
        let source = PathBuf::from("dante_data/analyses/2025-03-14-17-19-06_akjndkjenka_single".to_string());
        let motif_ids = vec!["ALS".to_string(), "DM2".to_string()];
        // dante_data/analyses/2025-03-14-17-19-06_akjndkjenka_single/vpuk-23-001440-A/data.json
        let data = Data { source, motif_ids };
        return ContentPage::SingleResults(data);
    }
}

fn view_header(mut content: Column<Message>, source: PathBuf) -> Column<Message> {
    content = content.push(row![
        container(button("Back").on_press(Message::Exit(source))).width(100),
        container(text("Result editor").size(App::H1_SIZE)).align_x(Horizontal::Center).width(Length::Fill),
        container(button("Save").on_press(Message::Save)).width(100).align_x(Horizontal::Right),
    ].padding(25).align_y(Vertical::Center));
    return content;
}

fn view_results<'a>(content: Column<'a, Message>, _data: &Data) -> Column<'a, Message> {
    let json = "dante_data/analyses/2025-03-14-17-19-06_akjndkjenka_single/vpuk-23-001440-A/data.json";
    let json: String = std::fs::read_to_string(json).expect("Cannot read file.");
    let json: Value = serde_json::from_str(&json).expect("JSON was not well-formatted");

    println!("{:?}", json);
    // let x = &json["x"];
    // for motif_id in &self.motif_ids {
    //     content = content.push(row![
    //         container(text(motif_id))
    //     ]);
    // }
    return content;
}
