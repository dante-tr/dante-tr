use iced::alignment::Horizontal;
use iced::widget::{button, column, horizontal_rule, scrollable};
use iced::{Element, Size};
use iced::widget::Column;
use iced::widget::{row, container, text};
use iced::Length;
use iced::alignment::Vertical;
use iced::widget::text_input;

use std::path::PathBuf;
use serde::{Serialize, Deserialize};
use std::io::{BufReader, BufRead};
use std::fs::File;
use std::iter::zip;
use std::path::Path;

use crate::ContentPage;
use crate::App;

#[derive(Debug, Clone)]
pub(crate) enum Message {
    Exit(PathBuf),
    Save,
    Edit(usize, String)
}

#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub(super) struct Data {
    source: PathBuf,
    meta_file: PathBuf,
    header: Vec<String>,
    content: Vec<String>,
    save_msg: String,
}

impl Data {
    pub(super) fn view(&self, _size: Size) -> Element<Message> {
        let mut content = column![].align_x(Horizontal::Center);

        content = view_header(content, self.source.clone(), &self.save_msg);
        content = view_metadata(content, &self.header, &self.content);
        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return content.into();
    }

    pub(super) fn update(&mut self, m: Message) {
        self.save_msg = "".to_string();
        use std::io::Write;
        match m {
            Message::Exit(_) => { unreachable!("Implemented in App::update."); }
            Message::Save => {
                let mut f = File::create(&self.meta_file).unwrap();
                for (h, c) in zip(&self.header, &self.content) {
                    writeln!(f, "{h}\t{c}").expect("Error writing.");
                }
                self.save_msg = "Saved! ".to_string();
            }
            Message::Edit(idx, entry) => { self.content[idx] = entry; }
        }
    }

    pub(super) fn open(source: PathBuf, meta_file: PathBuf) -> ContentPage {
        if !meta_file.exists() {
            std::fs::copy(App::get_filename("assets/includes/template.meta.tsv"), &meta_file).unwrap();
        }

        let (header, content) = read_meta_file(&meta_file);
        let save_msg = "".to_string();

        let data = Data {
            source,
            meta_file,
            header,
            content,
            save_msg
        };
        return ContentPage::MetadataEditor(data);
    }
}

pub(crate) fn read_meta_file(meta_file: &Path) -> (Vec<String>, Vec<String>) {
    let lines = BufReader::new(File::open(meta_file).expect("Cannot open metadata file.")).lines();
    let mut header = Vec::new();
    let mut content = Vec::new();
    for line in lines {
        let line = line.expect("Unable to read lines.");
        // TODO: make next lines more efficient
        let tmp: Vec<String> = line.split("\t").map(|x| x.to_string()).collect();
        assert!(tmp.len() == 2);
        header.push(tmp[0].clone());
        content.push(tmp[1].clone());
    }

    return (header, content);
}

fn view_metadata<'a>(content: Column<'a, Message>, keys: &'a[String], values: &'a[String]) -> Column<'a, Message> {
    let mut col = column![];

    let mut last_key_section = "";
    for (idx, (key, value)) in zip(keys, values).enumerate() {
        let tmp: Vec<&str> = key[1..].split(" - ").collect();
        let section = tmp[0];
        let entry = tmp[1];
        if last_key_section != section {
            col = col.push(row![
                container(horizontal_rule(1)),
                container(text(section)),
                container(horizontal_rule(1))
            ].padding(10).align_y(Vertical::Center));
            last_key_section = section;
        }

        let msg = move |x| Message::Edit(idx, x);
        col = col.push(row![
            container(text(entry)).width(200).align_x(Horizontal::Right).padding(App::PAD1),
            container(text_input("add value", value).on_input(msg)).padding(App::PAD1),
            container("").width(10)
        ].padding(10).align_y(Vertical::Center));
    }

    return content.push(scrollable(col));
}

fn view_header<'a>(mut content: Column<'a, Message>, source: PathBuf, save_msg: &'a str) -> Column<'a, Message> {
    content = content.push(row![
        container(button("Back").on_press(Message::Exit(source))).width(100),
        container(text("Metadata editor").size(App::H1_SIZE)).align_x(Horizontal::Center).width(Length::Fill),
        container(text(save_msg)).width(100).align_x(Horizontal::Right),
        container(button("Save").on_press(Message::Save)).align_x(Horizontal::Right),
    ].padding(25).align_y(Vertical::Center));
    return content;
}
