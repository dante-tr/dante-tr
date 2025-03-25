use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, column, container, horizontal_rule, horizontal_space, row, text};
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
    json: PathBuf,
}

impl Data {
    pub(super) fn view(&self, _size: Size) -> Element<Message> {
        let mut content = column![];

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
        let motif_ids: Vec<String> = [
            "ALS", "DM1", "DM2", "SCA27B", "SCA4"
        ].iter().map(|x| x.to_string()).collect();
        let source = "dante_data/analyses/2025-03-18-10-48-22_analysis1_single";
        let json = "dante_data/analyses/2025-03-18-10-48-22_analysis1_single/vpuk-23-001504-A/data_v2.json";
        // BAM???
        // meta???

        let source = PathBuf::from(source.to_string());
        let json = PathBuf::from(json.to_string());
        let data = Data { source, motif_ids, json };
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

fn view_results<'a>(mut content: Column<'a, Message>, data: &Data) -> Column<'a, Message> {
    let json: String = std::fs::read_to_string(&data.json).expect("Cannot read file.");
    let json: Value = serde_json::from_str(&json).expect("JSON was not well-formatted");

    content = content.push(horizontal_rule(1));
    content = content.push(row![container(text(format!("motif {} module {}", "ALS", 0)))]);

    content = content.push(row![
        horizontal_space(),
        container(text("Predicted results")),
        horizontal_space()
    ]);

    use iced::widget::text::Wrapping;
    // use iced::widget::text::Shaping;
    content = content.push(row![
        container(text("Repeat numbershdkfjsnkcdnskdcnskdnskdcnsdncskjdcnskdjcnskdjsn")).max_width(20),
        container(text("Repeat numbershdkfjsnkcdnskdcnskdnskdcnsdncskjdcnskdjcnskdjsn").wrapping(Wrapping::None)).clip(true),
    ].clip(true));

    // content = content.push(row![
    //     container(text("Sample")).width(Length::FillPortion(2)),
    //     container(text("Allele 1")).width(Length::FillPortion(6)),
    //     container(text("Allele 2")).width(Length::FillPortion(6)),
    //     container(text("Overall")).width(Length::FillPortion(7)),
    // ]);

    // content = content.push(row![
    //     container(text("No.")).width(Length::FillPortion(1)).clip(true),
    //     container(text("ID")).width(Length::FillPortion(1)).clip(true),

    //     container(text("Repeat number")).clip(true).width(Length::FillPortion(1)),
    //     container(text("Confidence")).width(Length::FillPortion(1)).clip(true),
    //     container(text("Pathogenicity")).width(Length::FillPortion(1)).clip(true),
    //     container(text("Spanning reads")).width(Length::FillPortion(1)).clip(true),
    //     container(text("Indel errors")).width(Length::FillPortion(1)).clip(true),
    //     container(text("Mismatch error")).width(Length::FillPortion(1)).clip(true),

    //     container(text("Repeat number")).width(Length::FillPortion(1)),
    //     container(text("Confidence")).width(Length::FillPortion(1)),
    //     container(text("Pathogenicity")).width(Length::FillPortion(1)),
    //     container(text("Spanning reads")).width(Length::FillPortion(1)),
    //     container(text("Indel errors")).width(Length::FillPortion(1)),
    //     container(text("Mismatch error")).width(Length::FillPortion(1)),

    //     container(text("???")).width(Length::FillPortion(1)),
    //     container(text("???")).width(Length::FillPortion(1)),
    //     container(text("???")).width(Length::FillPortion(1)),
    //     container(text("???")).width(Length::FillPortion(1)),
    //     container(text("???")).width(Length::FillPortion(1)),
    //     container(text("???")).width(Length::FillPortion(1)),
    //     container(text("???")).width(Length::FillPortion(1)),
    // ]);
    // for motif in json["motifs"].as_array().unwrap() {
    //     let x = motif["motif_id"].as_str().unwrap().to_string();
    //     if data.motif_ids.contains(&x) {
    //         content = content.push(horizontal_rule(2));
    //         content = content.push(row![container(text(motif["motif_id"].to_string()))]);
    //         content = content.push(row![container(text(motif["motif_stats"].to_string()))]);
    //         content = content.push(row![container(text(motif["phased_seqs"].to_string()))]);
    //         content = content.push(row![container(text(motif["nomenclatures"].to_string()))]);
    //         // println!("{:?}", motif["modules"]);
    //     }
    // }
    return content;
}
