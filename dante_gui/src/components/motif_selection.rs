use serde::{Serialize, Deserialize};
use iced::alignment::{Horizontal, Vertical};
use iced::widget::Row;
use iced::widget::checkbox;
use iced::widget::pick_list;
use iced::widget::tooltip;
use iced::widget::{column, container, horizontal_space, row, text};
use iced::{Element, Length, Size, Padding};
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::path::Path;
use std::error::Error;
use std::fs::File;
use native_dialog::FileDialog;

use crate::App;

#[derive(Debug, Clone)]
pub(crate) enum Message {
    SetMotifs(MotifFile),
    MotifCheckbox(usize, bool),
    MotifGroupbox(usize, bool),
}

#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub(crate) struct MotifSelection {
    pub(crate) selected: Option<MotifFile>,
    pub(crate) selected_file: Option<PathBuf>,
    pub(crate) motifs: Vec<(bool, String, Vec<String>, String)>,  // (checked, id, groups, description)
    pub(crate) groups: Vec<(bool, String)>,
}

impl MotifSelection {
    pub(crate) fn update(&mut self, m: Message) {
        match m {
            Message::SetMotifs(motif_file) => { update_motif_selection(self, motif_file) },
            Message::MotifGroupbox(idx, checked) => { toggle_group(self, idx, checked) },
            Message::MotifCheckbox(idx, checked) => { self.motifs[idx].0 = checked },
        }
    }

    pub(crate) fn view(&self, size: Size) -> Element<'_, Message> {
        return motif_part(self, size);
    }
}

fn motif_part(data: &MotifSelection, size: Size) -> Element<'_, Message>{
    let mut motif_part = column![].padding(5);
    let r0 = make_motif_selection(data.selected, &data.selected_file);

    let available_width = size.width as usize - 5 - 160 - 5 - 5;

    let mut i = 0;
    let mut r1 = row![].padding(5).align_y(Vertical::Center);
    r1 = r1.push(container(text("Group filter: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1));
    r1 = r1.extend(make_group_row(&data.groups, available_width, &mut i));
    r1 = r1.push(horizontal_space());
    // let r = std::convert::Into::<Element<Message>>::into(r).explain(iced::Color::BLACK);

    let mut i = 0;
    let mut r2 = row![].padding(5).align_y(Vertical::Center);
    r2 = r2.push(container(text("Motif filter: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1));
    r2 = r2.extend(make_checkbox_row(&data.motifs, available_width, &mut i));
    r2 = r2.push(horizontal_space());
    // let r = std::convert::Into::<Element<Message>>::into(r).explain(iced::Color::BLACK);

    // BUG: if available_width is too small, i is never increased
    let mut v: Vec<Element<_>> = Vec::new();
    while i < data.motifs.len() {
        let mut r3 = row![].padding(5).align_y(Vertical::Center);
        r3 = r3.push(container(text("")).width(160).align_x(Horizontal::Right));
        r3 = r3.extend(make_checkbox_row(&data.motifs, available_width, &mut i));
        r3 = r3.push(horizontal_space());
        // let r = std::convert::Into::<Element<Message>>::into(r).explain(iced::Color::BLACK);
        v.push(r3.into());
    }
    motif_part = motif_part.push(r0);
    motif_part = motif_part.push(r1);
    motif_part = motif_part.push(r2);
    motif_part = motif_part.extend(v);
    // let motif_part = std::convert::Into::<Element<Message>>::into(motif_part).explain(iced::Color::BLACK);
    return motif_part.into();
}

fn make_checkbox_row<'a>(motifs: &'a[(bool, String, Vec<String>, String)], available_width: usize, i: &mut usize) -> Vec<Element<'a, Message>> {
    const PAD: Padding = Padding { bottom: 0.0, top: 0.0, right: 15.0, left: 0.0 };
    let spacing = 15 /* checkbox */ + 10 /* between checkbox and label */ + 15 /* right padding */;
    let letter_width = 11;

    let mut v = Vec::new();
    let mut cur_width = 0;
    while *i < (*motifs).len() && cur_width + spacing + motifs[*i].1.len() * letter_width < available_width {
        let (ref checked, ref id, _, ref name) = &motifs[*i];
        let ii = *i;
        let f = move |b| Message::MotifCheckbox(ii, b);
        v.push(container(
            tooltip(
                checkbox(id, *checked).on_toggle(f),
                container(text(name)).padding(5).style(container::rounded_box),
                tooltip::Position::FollowCursor,
            )
        ).padding(PAD).into());
        cur_width += spacing + id.len() * letter_width;
        *i += 1;
    }
    return v;
}

fn make_group_row<'a>(groups: &'a[(bool, String)], available_width: usize, i: &mut usize) -> Vec<Element<'a, Message>> {
    const PAD: Padding = Padding { bottom: 0.0, top: 0.0, right: 15.0, left: 0.0};
    let spacing = 15 /*checkbox*/ + 10 /*between checkbox and label*/ + 15 /*right padding*/;
    let letter_width = 11;

    let mut v = Vec::new();
    let mut cur_width = 0;
    while *i < (*groups).len() && cur_width + spacing + groups[*i].1.len() * letter_width < available_width {
        let (ref checked, ref id) = &groups[*i];
        let ii = *i;
        let f = move |b| Message::MotifGroupbox(ii, b);
        v.push(container(checkbox(id, *checked).on_toggle(f)).padding(PAD).into());
        cur_width += spacing + id.len() * letter_width;
        *i += 1;
    }
    return v;
}

fn make_motif_selection(selected: Option<MotifFile>, selected_file: &Option<PathBuf>) -> Row<'_, Message> {
    let motif_files = [MotifFile::STRSet_20250311, MotifFile::Custom];

    let content = match selected {
        None => row![
            container(text("Motifs: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
            container(pick_list(motif_files, selected, Message::SetMotifs).placeholder("type")),
            container(text("")).width(Length::Fill).align_x(Horizontal::Left)
        ].padding(5).align_y(Vertical::Center),
        Some(MotifFile::Custom) => {
            let x: String = selected_file.clone().unwrap().to_string_lossy().to_string();
            row![
                container(text("Motifs: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
                container(pick_list(motif_files, selected, Message::SetMotifs).placeholder("type")).padding(App::PAD1),
                container(text(x)).width(Length::Fill).align_x(Horizontal::Left)
            ].padding(5).align_y(Vertical::Center)
        }
        Some(_) => {
            row![
                container(text("Motifs: ")).width(160).align_x(Horizontal::Right).padding(App::PAD1),
                container(pick_list(motif_files, selected, Message::SetMotifs).placeholder("type")),
                container(text("")).width(Length::Fill).align_x(Horizontal::Left)
            ].padding(5).align_y(Vertical::Center)
        }
    };
    return content;
}

fn update_motif_selection(data: &mut MotifSelection, motif_file: MotifFile) {
    use crate::analysis_common::{parse_motifs, get_groups};
    match motif_file {
        MotifFile::Custom => {
            if let Ok(Some(path)) = FileDialog::new().show_open_single_file() {
                data.selected = Some(motif_file);
                data.selected_file = Some(path);

                let format = validate_STR_format(data.selected_file.as_ref().unwrap());
                if format.is_ok() {
                    data.motifs = parse_motifs(data.selected_file.as_ref().unwrap());
                    data.groups = get_groups(data.motifs.as_ref());
                } else {
                    data.motifs = Vec::new();
                    data.groups = Vec::new();
                }
            }
        },
        _ => {
            let motif_str = motif_file.to_string();
            let path = PathBuf::from(App::DATA_DIR.to_string() + "/includes/" + &motif_str + ".tsv");
            data.selected = Some(motif_file);
            data.selected_file = Some(path);
            data.motifs = parse_motifs(data.selected_file.as_ref().unwrap());
            data.groups = get_groups(data.motifs.as_ref());
        }
    }
}

#[allow(non_snake_case)]
fn validate_STR_format(path: &Path) -> Result<(), Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let header = reader.lines().next().ok_or("File does not contain any lines.")??;
    let header: Vec<_> = header.split('\t').collect();

    if header[0] != "Disease ID"
        { return Err("1st column has incorrect name".into()); }
    if header[1] != "HGVS nomenclature (GRCh38 reference)"
        { return Err("2nd column has incorrect name".into()); }
    if header[2] != "Left flank"
        { return Err("3rd column has incorrect name".into()); }
    if header[3] != "Right flank"
        { return Err("4th column has incorrect name".into()); }
    if header[4] != "Groups"
        { return Err("5th column has incorrect name".into()); }
    if header[5] != "Disease name"
        { return Err("6th column has incorrect name".into()); }

    return Ok(());
}

fn toggle_group(data: &mut MotifSelection, idx: usize, checked: bool) {
    data.groups[idx].0 = checked;
    let group = data.groups[idx].1.clone();
    for x in &mut data.motifs { if x.2.contains(&group) { x.0 = checked; } }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) enum MotifFile {
    #[allow(non_camel_case_types)]
    STRSet_20220902,
    #[allow(non_camel_case_types)]
    STRSet_20250311,
    Custom,
}

impl std::fmt::Display for MotifFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::STRSet_20220902 => "STRSet_20220902",
            Self::STRSet_20250311 => "STRSet_20250311",
            Self::Custom => "custom",
        })
    }
}
