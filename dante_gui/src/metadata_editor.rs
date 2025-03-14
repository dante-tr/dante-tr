use iced::alignment::Horizontal;
use iced::widget::{button, column, scrollable};
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
    content: Vec<String>
}

impl Data {
    pub(super) fn view(&self, _size: Size) -> Element<Message> {
        let mut content = column![].align_x(Horizontal::Center);

        content = view_header(content, self.source.clone());
        content = view_metadata(content, &self.header, &self.content);
        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return content.into();
    }

    pub(super) fn update(&mut self, m: Message) {
        use std::io::Write;
        match m {
            Message::Exit(_) => { unreachable!("Implemented in App::update."); }
            Message::Save => {
                let mut f = File::create(&self.meta_file).unwrap();
                writeln!(f, "{}", self.header.join("\t")).expect("Error writing.");
                writeln!(f, "{}", self.content.join("\t")).expect("Error writing.");

            }
            Message::Edit(idx, entry) => { self.content[idx] = entry; }
        }
    }

    pub(super) fn open(source: PathBuf, meta_file: PathBuf) -> ContentPage {
        if !meta_file.exists() {
            std::fs::copy("./assets/template.meta.tsv", &meta_file).unwrap();
        }

        let mut lines = BufReader::new(File::open(&meta_file).expect("Cannot open metadata file.")).lines();
        let header = lines.next().unwrap().unwrap();
        let header: Vec<String> = header.split("\t").map(|x| x.to_string()).collect();
        let content = lines.next().unwrap().unwrap();
        let content: Vec<String> = content.split("\t").map(|x| x.to_string()).collect();

        let data = Data {
            source,
            meta_file,
            header,
            content
        };
        return ContentPage::MetadataEditor(data);
    }
}

fn view_metadata<'a>(content: Column<'a, Message>, keys: &'a[String], values: &'a[String]) -> Column<'a, Message> {
    let mut col = column![];

    for (idx, (key, value)) in zip(keys, values).enumerate() {
        let msg = move |x| Message::Edit(idx, x);
        col = col.push(row![
            container(text(key.clone())).width(350).align_x(Horizontal::Right).padding(App::PAD1),
            container(text_input("add value", value).on_input(msg)).padding(App::PAD1),
            container("").width(10)
        ].padding(10).align_y(Vertical::Center));
    }

    return content.push(scrollable(col));
}

fn view_header(mut content: Column<Message>, source: PathBuf) -> Column<Message> {
    content = content.push(row![
        container(button("Back").on_press(Message::Exit(source))).width(100),
        container(text("Metadata editor").size(App::H1_SIZE)).align_x(Horizontal::Center).width(Length::Fill),
        container(button("Save").on_press(Message::Save)).width(100).align_x(Horizontal::Right),
    ].padding(25).align_y(Vertical::Center));
    return content;
}
