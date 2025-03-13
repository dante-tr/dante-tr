use iced::alignment::Horizontal;
use iced::widget::{button, column, scrollable};
use iced::{Element, Size};

use std::path::PathBuf;
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::Write;

use crate::ContentPage;

#[derive(Debug, Clone)]
pub(crate) enum Message {
    Edit,
    Exit(PathBuf),
    SaveExit(PathBuf),
}

#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub(super) struct Data {
    path: PathBuf,
    bam_file: PathBuf,
}

impl Data {
    pub(super) fn view(&self, _size: Size) -> Element<Message> {
        let mut content = column![].align_x(Horizontal::Center);

        content = content.push(
            button("Exit").on_press(Message::Exit(self.path.clone()))
        );
        content = content.push(
            button("Save&Exit").on_press(Message::SaveExit(self.path.clone()))
        );
        // let content = std::convert::Into::<Element<Message>>::into(content).explain(iced::Color::BLACK);
        return scrollable(content).into();
    }

    pub(super) fn update(&mut self, m: Message) {
        match m {
            Message::Edit => {}
            Message::Exit(_) => { unreachable!() }
            Message::SaveExit(_) => { unreachable!() }
        }
    }

    // pub(super) fn init(path: PathBuf, analysis_name: String) -> ContentPage {
    //     let data = Data { path, analysis_name, ..Default::default() };
    //     data.save();
    //     // ContentPage::AnalysisSingle(data)
    //     ContentPage::default()
    // }

    fn save(&self) -> PathBuf {
        let json = serde_json::to_string(self).unwrap();
        let mut output = self.path.clone();
        output.push("params.json");
        let mut out = File::create(&output)
            .expect("Cannot open file for writing.");
        out.write_all(json.as_bytes())
            .expect("Cannot write to output file.");
        return output;
    }

    pub(super) fn load(mut path: PathBuf) -> Self {
        path.push("params.json");
        let json: String = std::fs::read_to_string(path)
            .expect("Cannot read file.");
        serde_json::from_str(&json)
            .expect("Cannot parse json.")
    }
}

pub(super) fn open(source: PathBuf, bam_file: PathBuf) -> ContentPage {
    println!("Edit metadata {:?} {:?}", source, bam_file);
    ContentPage::MetadataEditor(Data { path: source, bam_file })
}
