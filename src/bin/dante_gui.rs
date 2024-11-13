use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, column, row, text, image, text_input, container, checkbox, horizontal_rule};
use iced::{Element, Padding, Theme};
use native_dialog::FileDialog;
use std::path::PathBuf;

pub fn main() -> iced::Result {
    let settings = iced::window::Settings {
        size: iced::Size{width: 720.0, height: 480.0},
        ..Default::default()
    };
    iced::application("Dante", State::update, State::view)
        .window(settings)
        .theme(|_| { Theme::SolarizedDark })
        .run()
}

#[derive(Default)]
struct State {
    ref_file: Option<PathBuf>,
    bam_file: Option<PathBuf>,
    motif_file: Option<PathBuf>,
    output: Option<PathBuf>,
    out_bam: bool,
    correction: bool,
    dedup: bool,
    flank: usize,
    q: u8,
    score: Option<char>,
    print_quality: bool,
}

#[derive(Debug, Clone, Copy)]
enum Message {
    SelectFile,
    RunDante,
    Increment,
    Decrement,
    CheckboxToggled(bool),
}

impl State {
    fn update(&mut self, message: Message) {
        match message {
            Message::Increment => { },
            Message::Decrement => { },
            Message::SelectFile => { load_file(&mut self.bam_file); },
            Message::RunDante => {
                self.output = Some("kjsdkjn".into());
                opener::open("https://www.rust-lang.org").unwrap();
            },
            Message::CheckboxToggled(is_checked) => { self.out_bam = is_checked },
        }
    }
}

impl State {
    fn view(&self) -> Element<Message> {
        let tmp: &str = match self.bam_file.as_ref() {
            Some(x) => x.to_str().unwrap(),
            None => "Load some file"
        };
        let pad1 = Padding { right: 5.0, ..Padding::default() };
        let pad2 = Padding { left: 5.0, ..Padding::default() };
        let left_width = 120;
    
        column![
        column![
            image("assets/logo_cut.png").height(100),
        ].width(720.0).align_x(Horizontal::Right),
        horizontal_rule(2),
        column![
            row![
                container(text("Reference file:").line_height(2.2).size(12).width(left_width).align_x(Horizontal::Right)).padding(pad1),
                text_input("Type something here...", "tmp").on_input(tmp2).size(12),
                container(button("Load file").on_press(Message::SelectFile)).padding(pad2),
            ].padding(10.0).align_y(Vertical::Center),
            row![
                container(text("BAM file:").line_height(2.2).size(12).width(left_width).align_x(Horizontal::Right)).padding(pad1),
                text_input("Type something here...", "tmp").on_input(tmp2).size(12),
                container(button("Load file").on_press(Message::SelectFile)).padding(pad2),
            ].padding(10.0).align_y(Vertical::Center),
            row![
                container(text("Motif file:").line_height(2.2).size(12).width(left_width).align_x(Horizontal::Right)).padding(pad1),
                text_input("Type something here...", "tmp").on_input(tmp2).size(12),
                container(button("Load file").on_press(Message::SelectFile)).padding(pad2),
            ].padding(10.0).align_y(Vertical::Center),
            row![
                container(text("Output directory:").line_height(2.2).size(12).width(left_width).align_x(Horizontal::Right)).padding(pad1),
                text_input("Type something here...", "tmp").size(12),
                container(button("Load file").on_press(Message::SelectFile)).padding(pad2),
            ].padding(10.0).align_y(Vertical::Center),
            row![
                container("").width(left_width).padding(pad1),
                checkbox("Output BAM", self.out_bam).on_toggle(Message::CheckboxToggled),
                checkbox("Correction", self.correction).on_toggle(Message::CheckboxToggled),
                checkbox("Dedup", self.dedup).on_toggle(Message::CheckboxToggled),
                checkbox("Print quality", self.print_quality).on_toggle(Message::CheckboxToggled),
            ].padding(10.0).align_y(Vertical::Center),
            row![
                container("").width(left_width).padding(pad1),
                button("Run").on_press(Message::RunDante),
                text(tmp).size(12),
            ].padding(10.0).align_y(Vertical::Center),
            row![
                container(text("Result: ").line_height(2.2).size(12).align_x(Horizontal::Right)),
                button("Go to result").style(|theme: &Theme, status| {
                    let style = button::Style {
                        background: None,
                        text_color: iced::Color { r: 1.0, g: 1.0, b: 1.0, a: 1.0 },
                        border: iced::Border { color: iced::Color::BLACK, width: 1.0, radius: iced::border::Radius::default() },
                        shadow: iced::Shadow::default()
                    };
                    return style;
                }).on_press(Message::RunDante)
            ].padding(10.0).align_y(Vertical::Center),
        ].width(720.0).align_x(Horizontal::Left)
        ].into()
    }
}

fn load_file(result: &mut Option<PathBuf>) {
    let path = FileDialog::new()
        .add_filter("PNG Image", &["png"])
        .add_filter("JPEG Image", &["jpg", "jpeg"])
        .show_open_single_file()
        .unwrap();

    let path = match path {
        Some(path) => path,
        None => return,
    };

    *result = Some(path);
}

fn tmp2(_: String) -> Message {
    Message::Increment
}
