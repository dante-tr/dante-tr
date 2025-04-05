use iced::alignment::{Horizontal, Vertical};
use iced::widget::{
    button, column, container, horizontal_rule, horizontal_space, image, pick_list, row, scrollable, text, text_input, tooltip, vertical_space
};
use iced::widget::{Column, Row, Tooltip};
use iced::{Element, Length, Padding, Size};

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

const PADLR25: Padding = Padding { left: 25.0, right: 25.0, top: 0.0, bottom: 0.0 };

#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub(super) struct Data {
    source: PathBuf,
    motif_ids: Vec<String>,
    json: PathBuf,
    plots: PathBuf,
}

impl Data {
    pub(super) fn view(&self, _size: Size) -> Element<Message> {
        let mut content = column![];

        content = view_header(content, self.source.clone());
        for motif_id in &self.motif_ids {
            content = content.push(view_results(motif_id, self));
        }
        // for motif_id in motif_ids {}
        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return scrollable(content).into();
    }

    pub(super) fn update(&mut self, m: Message) {
        match m {
            Message::Exit(_) => { unreachable!("Implemented in App::update."); }
            Message::Save => { println!("Save"); }
        }
    }

    pub(super) fn open(motif_ids: Vec<String>, source: PathBuf, sample: String) -> ContentPage {
        let mut json = source.clone(); json.push(&sample); json.push("data_v2.json");
        let mut plots: PathBuf = source.clone(); plots.push(&sample); plots.push("plots");

        let data = Data { source, motif_ids, json, plots };
        return ContentPage::SingleResults(data);
    }
}

fn view_header(mut content: Column<Message>, source: PathBuf) -> Column<Message> {
    content = content.push(row![
        container(button("Back").on_press(Message::Exit(source))).width(100),
        container(text("Result editor").size(App::H1_SIZE)).align_x(Horizontal::Center).width(Length::Fill),
        container("").width(100).align_x(Horizontal::Right),
    ].padding(25).align_y(Vertical::Center));
    return content;
}

fn view_results<'a>(motif_id: &str, data: &Data) -> Column<'a, Message> {
    let mut motif_section = column![].padding(PADLR25);

    motif_section = motif_section.push(horizontal_rule(0));
    motif_section = motif_section.push(vertical_space().height(25));
    motif_section = motif_section.push(row![
        container(text(motif_id.to_string()).size(App::H1_SIZE)).padding(PADLR25),
        horizontal_space(),
        container(button("Save").on_press(Message::Save)).width(100).align_x(Horizontal::Right),
        horizontal_space().width(25)
    ]);
    motif_section = motif_section.push(vertical_space().height(25));

    motif_section = motif_section.push(view_predicted_table(motif_id, data));
    motif_section = motif_section.push(vertical_space().height(25));

    motif_section = motif_section.push(view_revised_table(motif_id, data));
    motif_section = motif_section.push(vertical_space().height(25));

    let plot_dir = data.plots.to_string_lossy().to_string();  // This is not correct.
    motif_section = motif_section.push(row![
        horizontal_space(),
        image(format!("{plot_dir}/{motif_id}_histogram.png")),
        horizontal_space()
    ]);

    motif_section = motif_section.push(vertical_space().height(25));
    return motif_section;
}

fn view_predicted_table<'a>(motif_id: &str, data: &Data) -> Column<'a, Message> {
    let mut predicted_table = column![].padding(PADLR25);
    predicted_table = predicted_table.extend(view_table_header("Predicted results".to_string()));
    let json: String = std::fs::read_to_string(&data.json).expect("Cannot read file.");
    let json: Value = serde_json::from_str(&json).expect("JSON was not well-formatted");

    for motif in json["motifs"].as_array().unwrap() {
        let motif_id2 = motif["motif_id"].as_str().unwrap().to_string();
        if motif_id2 == motif_id {
            for (i, module) in motif["modules"].as_array().unwrap().iter().enumerate() {
                predicted_table = predicted_table.push(horizontal_rule(0));
                predicted_table = predicted_table.push(view_motif_module(motif, module, i));
            }
        }
    }
    predicted_table = predicted_table.push(horizontal_rule(0));
    return predicted_table;
}

fn view_revised_table<'a>(motif_id: &str, data: &Data) -> Column<'a, Message> {
    let mut predicted_table = column![].padding(PADLR25);
    predicted_table = predicted_table.extend(view_table_header("Revised results".to_string()));
    let json: String = std::fs::read_to_string(&data.json).expect("Cannot read file.");
    let json: Value = serde_json::from_str(&json).expect("JSON was not well-formatted");

    for motif in json["motifs"].as_array().unwrap() {
        let motif_id2 = motif["motif_id"].as_str().unwrap().to_string();
        if motif_id2 == motif_id {
            for (i, module) in motif["modules"].as_array().unwrap().iter().enumerate() {
                predicted_table = predicted_table.push(horizontal_rule(0));
                predicted_table = predicted_table.push({
                    let module_idx = i;
                    let a1_row = view_revised_allele_row(&module["allele_1"], 1);
                    let a2_row = view_revised_allele_row(&module["allele_2"], 2);
                    let stats_row = view_stats_row(module);

                    let mut content = column![];
                    content = content.push(row![
                        container(text(module_idx.to_string())).width(Length::FillPortion(1)),
                        column![
                            a1_row,
                            a2_row
                        ].width(Length::FillPortion(7)),
                        stats_row
                    ].align_y(Vertical::Center));
                    content
                });
            }
        }
    }
    predicted_table = predicted_table.push(horizontal_rule(0));
    return predicted_table;
}

fn view_table_header<'a>(title: String) -> Vec<Element<'a, Message>> {
    return vec![
        row![
            container(text("")).width(Length::FillPortion(1)),
            container(text(title)).width(Length::FillPortion(7)),
            container(text("Overall")).width(Length::FillPortion(7)),
        ].into(),
        row![
            view_header_cell("Mod", "Module"),

            view_header_cell("A", "Allele"),
            view_header_cell("Pred", "Prediction"),
            view_header_cell("Conf", "Confidence"),
            view_header_cell("Pat", "Pathogenicity"),
            view_header_cell("SR", "Spanning reads"),
            view_header_cell("I/D", "Indel errors"),
            view_header_cell("X", "Mismatch errors"),

            view_header_cell("Conf", "Confidence"),
            view_header_cell("BDR", "BAM-Dante reads"),
            view_header_cell("RT", "Reads total"),
            view_header_cell("RS", "Reads spanning"),
            view_header_cell("RP", "Reads partial"),
            view_header_cell("I/D", "Indel errors"),
            view_header_cell("X", "Mismatch errors"),
        ].into()
    ];
}

fn view_header_cell<'a>(short: &'a str, long: &'a str) -> Tooltip<'a, Message> {
    tooltip(
        // TODO: how to clip text instead of creating second line?
        container(text(short)).width(Length::FillPortion(1)),
        container(text(long)).padding(5).style(container::rounded_box),
        tooltip::Position::FollowCursor
    )
}

fn view_allele_row<'a>(allele_data: &Value, num: usize) -> Row<'a, Message> {
    let pat: &[_] = &['"', ' '];
    let pred = allele_data[0].to_string().trim_matches(pat).to_string();
    let conf = allele_data[1].to_string().trim_matches(pat).to_string();
    let indels = allele_data[2].to_string().trim_matches(pat).to_string();
    let matches = allele_data[3].to_string().trim_matches(pat).to_string();

    row![
        container(text(num.to_string()))           .width(Length::FillPortion(1)),
        container(text(pred)).width(Length::FillPortion(1)),
        container(text(conf)).width(Length::FillPortion(1)),
        container(text("Benign"))                  .width(Length::FillPortion(1)),
        container(text(allele_data[4].to_string())).width(Length::FillPortion(1)),
        container(text(indels)).width(Length::FillPortion(1)),
        container(text(matches)).width(Length::FillPortion(1)),
    ]
        .padding(Padding {left: 0.0, top: 10.0, right: 0.0, bottom: 10.0})
        .align_y(Vertical::Center)
}

fn view_revised_allele_row<'a>(allele_data: &Value, num: usize) -> Row<'a, Message> {
    let pat: &[_] = &['"', ' '];
    let pred = allele_data[0].to_string().trim_matches(pat).to_string();
    let conf = allele_data[1].to_string().trim_matches(pat).to_string();
    let indels = allele_data[2].to_string().trim_matches(pat).to_string();
    let matches = allele_data[3].to_string().trim_matches(pat).to_string();

    let value = "E";
    let selected: Option<Status> = None;
    let options = [ Status::Benign, Status::Affected ];
    let on_selected = |_| { Message::Save };


    row![
        container(text(num.to_string()))           .width(Length::FillPortion(1)),
        container(text_input(&pred, value).width(30))        .width(Length::FillPortion(1)),
        container(text(conf))                      .width(Length::FillPortion(1)),
        container(pick_list(options, selected, on_selected)).width(Length::FillPortion(1)),
        container(text(allele_data[4].to_string())).width(Length::FillPortion(1)),
        container(text(indels))                    .width(Length::FillPortion(1)),
        container(text(matches))                   .width(Length::FillPortion(1)),
    ]
        .padding(Padding {left: 0.0, top: 5.0, right: 0.0, bottom: 5.0})
        .align_y(Vertical::Center)
}

fn view_stats_row<'a>(module: &Value) -> Row<'a, Message> {
    let pat: &[_] = &['"', ' '];
    let conf = module["stats"][0].to_string().trim_matches(pat).to_string();
    let indels = module["stats"][1].to_string().trim_matches(pat).to_string();
    let mismatches = module["stats"][2].to_string().trim_matches(pat).to_string();
    let spanning = module["reads_spanning"].to_string();
    let flanking = module["reads_flanking"].to_string();
    row![
        container(text(conf)).width(Length::FillPortion(1)),
        container(text("???")).width(Length::FillPortion(1)),
        container(text("???")).width(Length::FillPortion(1)),
        container(text(spanning)).width(Length::FillPortion(1)),
        container(text(flanking)).width(Length::FillPortion(1)),
        container(text(indels)).width(Length::FillPortion(1)),
        container(text(mismatches)).width(Length::FillPortion(1)),
    ].width(Length::FillPortion(7))
}

fn view_motif_module<'a, 'b>(
    motif: &'b Value, module: &'b Value, module_idx: usize
) -> Column<'a, Message> {
    let a1_row = view_allele_row(&module["allele_1"], 1);
    let a2_row = view_allele_row(&module["allele_2"], 2);
    let stats_row = view_stats_row(module);

    let mut content = column![];
    content = content.push(row![
        container(text(module_idx.to_string())).width(Length::FillPortion(1)),
        column![
            a1_row,
            a2_row
        ].width(Length::FillPortion(7)),
        stats_row
    ].align_y(Vertical::Center));
    return content;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
enum Status {
    Benign,
    Affected
}

impl std::fmt::Display for Status {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Benign => "Benign",
            Self::Affected => "Affected",
        })
    }
}
