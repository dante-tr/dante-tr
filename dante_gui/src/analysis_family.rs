use iced::widget::column;
use iced::Element;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Data { 
    pub analysis_name: String
}

#[derive(Debug, Clone)]
pub(crate) enum Message { }


pub fn view(data: &Data) -> Element<Message> {
    println!("{:?}", data);
    column![].into()
}
