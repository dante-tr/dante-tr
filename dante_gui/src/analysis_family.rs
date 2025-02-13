use iced::widget::column;
use iced::Element;

use crate::App;
use crate::Message;

pub fn view(_state: &App) -> Element<Message> {
    column![].into()
}
