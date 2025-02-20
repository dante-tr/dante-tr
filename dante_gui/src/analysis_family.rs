use std::path::PathBuf;

use iced::widget::column;
use iced::Element;

use crate::ContentPage;

#[derive(Debug, Clone)]
pub(crate) enum Message {
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Data { 
    pub path: PathBuf,
    pub analysis_name: String
}

impl Data {
    pub(crate) fn init(path: PathBuf, analysis_name: String) -> ContentPage {
        ContentPage::AnalysisFamily(Data { path, analysis_name })
    }

    pub(crate) fn view(&self) -> Element<Message> {
        println!("{:?}", self);
        column![].into()
    }

    pub(crate) fn update(&mut self, m: Message) {
        println!("{:?}\n{:?}", self, m)
    }
}
