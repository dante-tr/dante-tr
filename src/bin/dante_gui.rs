use std::path::PathBuf;
use iced::Element;
use iced::widget::{button, column, row, text, image};
use native_dialog::FileDialog;

pub fn main() -> iced::Result {
    let settings = iced::window::Settings {
        size: iced::Size{width: 800.0, height: 600.0},
        ..Default::default()
    };
    iced::application("Dante", State::update, State::view)
        .window(settings)
        .run()
}

#[derive(Debug, Clone, Copy)]
enum Message {
    SelectFile,
    RunDante,
    Increment,
    Decrement,
}

#[derive(Default)]
struct State {
    bam_file: Option<PathBuf>,
    hgvs_nomenclature: Option<PathBuf>,
    result_report: Option<PathBuf>,
    value: i64,
}

impl State {
    fn update(&mut self, message: Message) {
        match message {
            Message::Increment => {
                self.value += 1;
            }
            Message::Decrement => {
                self.value -= 1;
            }
            Message::SelectFile => {
                let path = FileDialog::new()
                    // .set_location("~/Desktop")
                    .add_filter("PNG Image", &["png"])
                    .add_filter("JPEG Image", &["jpg", "jpeg"])
                    .show_open_single_file()
                    .unwrap();

                let path = match path {
                    Some(path) => path,
                    None => return,
                };

                self.bam_file = Some(path);
            },
            Message::RunDante => {
                self.result_report = Some("kjsdkjn".into());
            }
        }
    }

    fn view(&self) -> Element<Message> {
        let tmp: &str = match self.bam_file.as_ref() {
            Some(x) => x.to_str().unwrap(),
            None => "Load some file"
        };
        column![
            row![
                image("logo_cut.png").height(150)
            ],
            row![
                column![
                    text("skjbjh").size(12),
                ].padding(20),
                column![
                    button("Increment").on_press(Message::Increment),
                    text(self.value).size(12),
                    button("Decrement").on_press(Message::Decrement),
                    text(tmp).size(12),
                    button("Load file").on_press(Message::SelectFile),
                    button("Run").on_press(Message::RunDante),
                ]
                .padding(20),
            ]
        ].into()
    }
}
