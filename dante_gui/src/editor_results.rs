use iced::widget::{column, scrollable};
use iced::widget::text_editor;
use iced::{Element, Size};

use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;
use serde::{Serialize, Deserialize};
use serde_json::Value;

use crate::ContentPage;

#[derive(Debug, Clone)]
pub(super) enum Message {
    Exit(PathBuf),
    RevisionChanged(String, usize, usize, usize),   // text, motif_pos, module_idx, allele
    PatChanged(Status, usize, usize, usize),        // ^-
    Save(usize),
    SaveAll,

    ActionPerformed(text_editor::Action),

    // motif specific actions
    QCToggle(usize, bool),
    QCEdit(usize, String),
    InterpretationEdit(usize, String),
    ShowAlignments(String),
}

#[derive(Debug, Default)]
pub(super) struct Data {
    source: PathBuf,
    sample: String,
    json: PathBuf,
    plots: PathBuf,

    motif_ids: Vec<String>,
    motif_names: Vec<String>,
    revisions: Vec<Revision>,

    interpretation: text_editor::Content,
    save_msg: String,
}

impl Data {
    pub(super) fn view(&self, size: Size) -> Element<Message> {
        let mut content = column![];

        content = content.push(view::header(self.source.clone(), self));
        content = content.push(view::general_data(self));
        for motif_idx in 0..self.motif_ids.len() {
            content = content.push(view::motif_data(self, motif_idx, size));
        }
        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return scrollable(content).into();
    }

    pub(super) fn update(&mut self, m: Message) {
        self.save_msg = "".to_string();
        match m {
            Message::Exit(_) => { unreachable!("Implemented in App::update."); }
            Message::SaveAll => { self.save_all(); self.save_msg = "Saved! ".to_string(); }
            Message::Save(motif_idx) => { self.save_motif(motif_idx); }
            Message::PatChanged(status, motif_idx, module_idx, allele_idx) => {
                self.revisions[motif_idx].modules[module_idx][allele_idx].2 = Some(status);
            }
            Message::RevisionChanged(text, motif_idx, module_idx, allele_idx) => {
                self.revisions[motif_idx].modules[module_idx][allele_idx].0 = text;
            }
            Message::ShowAlignments(motif_id) => {
                let mut alignment = self.source.clone();
                alignment.push(&self.sample);
                alignment.push("alignments");
                alignment.push(format!("{}.html", motif_id));
                opener::open(alignment).unwrap();
            }
            Message::ActionPerformed(x) => {
                self.interpretation.perform(x);
            }
            Message::QCToggle(idx, value) => {
                self.revisions[idx].qc_passed = value;
            }
            Message::QCEdit(idx, value) => {
                self.revisions[idx].qc_notes = value;
            }
            Message::InterpretationEdit(idx, value) => {
                self.revisions[idx].locus_interpretation = value;
            }
        }
    }

    pub(super) fn open(
        motif_ids: Vec<String>, motif_names: Vec<String>, source: PathBuf, sample: String
    ) -> ContentPage {
        let mut json = source.clone(); json.push(&sample); json.push("data_v2.json");
        let mut plots: PathBuf = source.clone(); plots.push(&sample); plots.push("plots");

        let mut revisions: Vec<Revision> = Vec::with_capacity(motif_ids.len());
        for motif_id in &motif_ids {
            let mut revision = source.clone(); revision.push("revisions"); revision.push(format!("{}.json", motif_id));
            if revision.exists() {
                let rev_json: String = std::fs::read_to_string(revision).expect("Cannot read file.");
                let rev = serde_json::from_str(&rev_json).expect("Cannot parse json.");
                revisions.push(rev);
            } else {
                revisions.push(Revision::empty(&json, motif_id));
            }
        }

        let mut interpretation_path = source.clone();
        interpretation_path.push("revisions");
        interpretation_path.push("interpretation.txt");

        let interpretation = match interpretation_path.exists() {
            false => { text_editor::Content::new() }
            true => {
                let mut f = File::open(interpretation_path).expect("Cannot open interpretation");
                let mut buf = String::new();
                f.read_to_string(&mut buf).expect("Cannot read interpretation");
                text_editor::Content::with_text(&buf)
            },
        };
        let save_msg = "Saved! ".to_string();

        let data = Data {
            source,
            sample,
            json,
            plots,

            motif_ids,
            motif_names,
            revisions,

            interpretation,
            save_msg,
        };
        return ContentPage::SingleResults(data);
    }

    fn save_motif(&mut self, idx: usize) {
        let mut output = self.source.clone(); output.push("revisions");
        let _ = std::fs::create_dir(&output); /* ignoring error when dir already exists */
        output.push(format!("{}.json", self.motif_ids[idx]));

        let mut out = File::create(&output).expect("Cannot open file for writing.");
        self.revisions[idx].inherit_defaults();
        let json = serde_json::to_string(&self.revisions[idx]).unwrap();
        out.write_all(json.as_bytes()).expect("Cannot write to output file.");
    }

    fn save_all(&mut self) {
        for idx in 0..self.motif_ids.len() { self.save_motif(idx); }

        let mut output = self.source.clone(); output.push("revisions");
        let _ = std::fs::create_dir(&output); /* ignoring error when dir already exists */
        output.push("interpretation.txt");

        let mut out = File::create(&output).expect("Cannot open file for writing.");
        out.write_all(self.interpretation.text().as_bytes()).expect("Cannot write to output file.");
    }
}

mod view {
    use super::{Message, Data, Status};

    use iced::alignment::{Horizontal, Vertical};
    use iced::widget::{
        button, checkbox, column, container, horizontal_rule, horizontal_space,
        image, pick_list, row, text, text_editor, text_input, tooltip, vertical_space
    };
    use iced::widget::{Column, Row, Tooltip, Space};
    use iced::{Element, Length, Padding, Size};
    use std::path::PathBuf;
    use serde_json::Value;

    use crate::App;

    const PADLR25: Padding = Padding { top: 0.0, right: 25.0, bottom: 0.0, left: 25.0 };
    const PADLRB25: Padding = Padding { top: 0.0, right: 25.0, bottom: 25.0, left: 25.0 };
    const PADTB5: Padding = Padding { top: 5.0, right: 0.0, bottom: 5.0, left: 0.0 };
    const PADT5: Padding = Padding { top: 5.0, right: 0.0, bottom: 0.0, left: 0.0 };

    // this fn could be as follows, if I will have some problems with coupled lifetimes
    // pub(super) fn header<'a, 'b: 'a>(source: PathBuf, data: &'b Data) -> Row<'a, Message>
    pub(super) fn header(source: PathBuf, data: &Data) -> Row<Message> {
        row![
            container(button("Back").on_press(Message::Exit(source))).width(100),
            container(text("Result editor").size(App::H1_SIZE)).align_x(Horizontal::Center).width(Length::Fill),
            container(text(&data.save_msg)).align_x(Horizontal::Right).width(100),
            container(button("Save all").on_press(Message::SaveAll)).align_x(Horizontal::Right),
        ].padding(25).align_y(Vertical::Center)
    }

    pub(super) fn general_data(data: &Data) -> Element<Message> {
        let bam_line = format!("BAM ID: {}", data.sample);

        let mut general_data = column![].padding(PADLR25);
        general_data = general_data.push(horizontal_rule(0));
        general_data = general_data.push(vertical_space().height(25));
        general_data = general_data.push(container(text(bam_line).size(App::H1_SIZE)).padding(PADLRB25));

        general_data = general_data.push(container(text("Interpretation of results")).padding(PADLR25));
        general_data = general_data.push(container(
            text_editor(&data.interpretation).placeholder("Interpretation of results").on_action(Message::ActionPerformed)
        ).padding(PADLRB25));
        general_data.into()
    }

    pub(super) fn motif_data<'a>(data: &Data, motif_idx: usize, size: Size) -> Column<'a, Message> {
        let json: String = std::fs::read_to_string(&data.json).expect("Cannot read file.");
        let json: Value = serde_json::from_str(&json).expect("JSON was not well-formatted");
        let motif_id: &str = &data.motif_ids[motif_idx];
        let motif: &Value = json["motifs"].as_array().unwrap().iter().find(|x| x["motif_id"] == motif_id).unwrap();

        let mut motif_section = column![].padding(PADLR25);

        motif_section = motif_section.push(horizontal_rule(0));
        motif_section = motif_section.push(vertical_space().height(25));
        let motif_text = data.motif_ids[motif_idx].to_string() + " - " + &data.motif_names[motif_idx];
        motif_section = motif_section.push(row![
            container(text(motif_text).size(App::H1_SIZE)).padding(PADLR25),
            horizontal_space(),
            container(button("Save").on_press(Message::Save(motif_idx))).width(100).align_x(Horizontal::Right),
            horizontal_space().width(25)
        ]);
        let modules = view_modules(&motif["motif_stats"]["modules"]);
        motif_section = motif_section.push(row![
            container(text(format!("modules: {modules}"))).padding(PADLR25)
        ]);

        motif_section = motif_section.push(vertical_space().height(25));

        motif_section = motif_section.push(view_predicted_table(&data.motif_ids[motif_idx], data, motif));
        motif_section = motif_section.push(vertical_space().height(25));

        motif_section = motif_section.push(view_revised_table(&data.motif_ids[motif_idx], data, motif));
        motif_section = motif_section.push(vertical_space().height(25));

        {  // add_plots
            let plot_dir = data.plots.to_string_lossy().to_string();  // This is not correct.
            // let motif_id = &data.motif_ids[motif_idx];

            // let motif = json["motifs"].as_array().unwrap().iter().find(|x| x["motif_id"] == *motif_id).unwrap();

            for module_pos in 0..motif["modules"].as_array().unwrap().len() {
                motif_section = motif_section.push(row![
                    horizontal_space(),
                    image(format!("{plot_dir}/{motif_id}_{module_pos}_histogram.png")),
                    Space::with_width(50),
                    right_panel(data, motif_idx, &motif["modules"][module_pos], size),
                    horizontal_space()
                ].align_y(Vertical::Center));
                motif_section = motif_section.push(vertical_space().height(25));
            }
        }

        return motif_section;
    }

    fn right_panel<'a>(data: &Data, motif_idx: usize, module: &Value, size: Size) -> Column<'a, Message> {
        let pat: &[_] = &['"', ' '];
        let mut res: Vec<Element<'a, Message>> = Vec::new();
        res.push(text("Nomenclature counts (5 most common):").into());
        for nomenclature in module["nomenclatures"].as_array().unwrap() {
            res.push(text(format!("{}x {}",
                nomenclature["count"],
                nomenclature["noms"][0].to_string().trim_matches(pat).replace("]", "] ")
            )).into());
        }

        let tmp = &data.revisions[motif_idx];
        let toggle_message = move |x| { Message::QCToggle(motif_idx, x) };
        res.push(container(checkbox("QC passed", tmp.qc_passed).on_toggle(toggle_message)).padding(PADTB5).into());
        let input_message = move |x| { Message::QCEdit(motif_idx, x) };
        res.push(container(text("QC reason")).padding(PADT5).into());
        res.push(container(text_input("QC reason", &tmp.qc_notes).on_input(input_message)).padding(PADTB5).into());
        let input_message2 = move |x| { Message::InterpretationEdit(motif_idx, x) };
        res.push(container(text("locus interpretation")).padding(PADT5).into());
        res.push(container(text_input("locus interpretation", &tmp.locus_interpretation).on_input(input_message2)).padding(PADTB5).into());

        let m = Message::ShowAlignments(data.motif_ids[motif_idx].to_string());
        res.push(container(button("View alignments").on_press(m)).padding(PADTB5).into());
        let w = size.width - 640.0 /* plot */ - 4.0 * 25.0 /* margins */ - 50.0 /* mid sep */;
        return column![].extend(res).width(w).max_width(800);
    }

    fn view_modules(modules: &Value) -> String {
        let mods = modules.as_array().unwrap();
        let result = mods.iter()
            .skip(1).take(mods.len() - 2)  // skip flanks (first and last)
            .map(|m| format!("{}[{}]", m[0].as_str().unwrap(), m[1].as_u64().unwrap()))
            .collect::<Vec<String>>().join("  ");
        return result;
    }

    fn view_predicted_table<'a>(_motif_id: &str, _data: &Data, motif: &Value) -> Column<'a, Message> {
        // let json: String = std::fs::read_to_string(&data.json).expect("Cannot read file.");
        // let json: &Value = &serde_json::from_str(&json).expect("JSON was not well-formatted");
        // let motif_pos = data.motif_ids.iter().position(|x| x == motif_id).unwrap();
        // let motif: &Value = json["motifs"].as_array().unwrap().iter().find(|x| x["motif_id"] == motif_id).unwrap();

        let mut predicted_table = column![].padding(PADLR25);
        predicted_table = predicted_table.extend(view_table_header("Predicted results".to_string()));

        for (i, module) in motif["modules"].as_array().unwrap().iter().enumerate() {
            predicted_table = predicted_table.push(horizontal_rule(0));
            predicted_table = predicted_table.push(view_motif_module(module, i));
        }
        predicted_table = predicted_table.push(horizontal_rule(0));
        return predicted_table;
    }

    fn view_revised_table<'a>(motif_id: &str, data: &Data, motif: &Value) -> Column<'a, Message> {
        // let json: String = std::fs::read_to_string(&data.json).expect("Cannot read file.");
        // let json: &Value = &serde_json::from_str(&json).expect("JSON was not well-formatted");
        // let motif: &Value = json["motifs"].as_array().unwrap().iter().find(|x| x["motif_id"] == motif_id).unwrap();

        let motif_pos = data.motif_ids.iter().position(|x| x == motif_id).unwrap();

        let mut predicted_table = column![].padding(PADLR25);
        predicted_table = predicted_table.extend(view_table_header("Revised results".to_string()));

        for module_pos in 0..motif["modules"].as_array().unwrap().len() {
            predicted_table = predicted_table.push(horizontal_rule(0));
            predicted_table = predicted_table.push(view_motif_module_revised(data, motif, motif_pos, module_pos));
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
            container(text("unk"))                  .width(Length::FillPortion(1)),
            container(text(allele_data[4].to_string())).width(Length::FillPortion(1)),
            container(text(indels)).width(Length::FillPortion(1)),
            container(text(matches)).width(Length::FillPortion(1)),
        ]
            .padding(Padding {left: 0.0, top: 10.0, right: 0.0, bottom: 10.0})
            .align_y(Vertical::Center)
    }

    fn view_revised_allele_row<'a>(
        data: &Data, allele_data: &Value, motif_pos: usize, module_idx: usize, num: usize
    ) -> Row<'a, Message> {
        let pat: &[_] = &['"', ' '];
        let conf    = allele_data[1].to_string().trim_matches(pat).to_string();
        let indels  = allele_data[2].to_string().trim_matches(pat).to_string();
        let matches = allele_data[3].to_string().trim_matches(pat).to_string();

        let pred = allele_data[0].to_string().trim_matches(pat).to_string();
        let value = &data.revisions[motif_pos].modules[module_idx][num-1].0;
        let on_input_f = move |x| { Message::RevisionChanged(x, motif_pos, module_idx, num-1) };
        let txt_in = text_input(&pred, value).on_input(on_input_f).width(40);

        let options = [
            Status::Benign,
            Status::LikelyBenign,
            Status::Premutation,
            Status::LikelyPathogenic,
            Status::Pathogenic,
            Status::Unknown
        ];
        let selected: Option<Status> = data.revisions[motif_pos].modules[module_idx][num-1].2;
        let on_pick_f = move |x| { Message::PatChanged(x, motif_pos, module_idx, num-1) };
        let pck_lst = pick_list(options, selected, on_pick_f);

        // type VElem<'a> = Vec<Element<'a, Message>>;
        let x: Vec<Element<'a, Message>> = vec![
            text(num.to_string()).into(),
            txt_in.into(),
            text(conf).into(),
            pck_lst.into(),
            text(allele_data[4].to_string()).into(),
            text(indels).into(),
            text(matches).into()
        ];

        let result = row![]
            .padding(Padding {left: 0.0, top: 5.0, right: 0.0, bottom: 5.0})
            .align_y(Vertical::Center);

        let allele_row: Vec<Element<'a, Message>> = x.into_iter()
            .map(|y| container(y).width(Length::FillPortion(1)).into()).collect();
        result.extend(allele_row)
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

    fn view_motif_module_revised<'a>(data: &Data, motif: &Value, motif_pos: usize, module_idx: usize) -> Column<'a, Message> {
        let module = &motif["modules"][module_idx];

        let a1_row = view_revised_allele_row(data, &module["allele_1"], motif_pos, module_idx, 1);
        let a2_row = view_revised_allele_row(data, &module["allele_2"], motif_pos, module_idx, 2);
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
    }

    fn view_motif_module<'a>(
        module: &Value, module_idx: usize
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
}

mod update {}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub(super) enum Status {
    Benign,
    Premutation,
    Pathogenic,
    #[default]
    Unknown,
    LikelyBenign,
    LikelyPathogenic,
}

impl std::fmt::Display for Status {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Benign => "ben",
            Self::Premutation => "pre",
            Self::Pathogenic => "pat",
            Self::Unknown => "unk",
            Self::LikelyPathogenic => "lpat",
            Self::LikelyBenign => "lben",
        })
    }
}

impl Status {
    pub(super) fn to_typst(self) -> String {
        match self {
            Self::Benign => "benign".to_string(),
            Self::Premutation => "premutation".to_string(),
            Self::Pathogenic => "pathogenic".to_string(),
            Self::Unknown => "unknown".to_string(),
            Self::LikelyPathogenic => "likely_pathogenic".to_string(),
            Self::LikelyBenign => "likely_benign".to_string()
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub(crate) struct Revision {
    pub(crate) motif_id: String,
    pub(crate) modules: Vec<[(String /* rev_pred */, String /* raw_predict */, Option<Status> /* pathogenicity */); 2]>,
    pub(crate) qc_passed: bool,
    pub(crate) qc_notes: String,
    pub(crate) locus_interpretation: String
}

impl Revision {
    fn empty(json_path: &PathBuf, motif_id: &str) -> Self {
        let json: String = std::fs::read_to_string(json_path).expect("Cannot read file.");
        let json: Value = serde_json::from_str(&json).expect("Cannot parse json.");
        let motif = json["motifs"].as_array().unwrap().iter().find(|x| x["motif_id"] == motif_id).unwrap();
        let modules = motif["modules"].as_array().unwrap();

        let mut v = Vec::new();

        let get_allele = |x: &Value| {
            if x.is_string() { x.as_str().unwrap().to_owned() }
            else { x.to_string() }
        };

        for m in modules {
            v.push([
                ("".to_string(), get_allele(&m["allele_1"][0]), Some(Status::Unknown)),
                ("".to_string(), get_allele(&m["allele_2"][0]), Some(Status::Unknown))
            ]);
        }

        return Self {
            motif_id: motif_id.to_string(),
            modules: v,
            ..Default::default()
        };
    }

    fn inherit_defaults(&mut self) {
        for m in &mut self.modules {
            if m[0].0.is_empty() { m[0].0 = m[0].1.clone() }
            if m[1].0.is_empty() { m[1].0 = m[1].1.clone() }
        }
    }
}
