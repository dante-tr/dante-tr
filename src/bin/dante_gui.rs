use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, checkbox, column, container, horizontal_rule, image, row, text, text_input, Row};
use iced::{Element, Padding, Theme};
use native_dialog::FileDialog;
use std::path::PathBuf;
use remastr::run;
use std::process::Command;

pub fn main() -> iced::Result {
    let settings = iced::window::Settings {
        size: iced::Size{width: 720.0, height: 480.0},
        ..Default::default()
    };
    iced::application("Dante", State::update, State::view)
        .window(settings)
        .theme(|_| { Theme::Nord })
        .run()
}

#[derive(Debug, Default)]
struct State {
    ref_file: Option<PathBuf>,
    bam_file: Option<PathBuf>,
    motif_file: Option<PathBuf>,
    output: Option<PathBuf>,
    out_bam: bool,
    correction: bool,
    dedup: bool,
    _flank: usize,
    _q: u8,
    _score: Option<char>,
    print_quality: bool,
}

#[derive(Debug, Clone)]
enum Message {
    RefChanged(String),
    SelectRef,
    BamChanged(String),
    SelectBam,
    MotifChanged(String),
    SelectMotif,
    OutdirChanged(String),
    SelectOutdir,

    RunDante,
    OpenResults,
    CheckboxToggled(bool),
}

impl State {
    fn update(&mut self, message: Message) {
        match message {
            Message::RefChanged(content) => { self.ref_file = Some(PathBuf::from(content)); },
            Message::SelectRef => { Self::load_file(&mut self.ref_file); },
            Message::BamChanged(content) => { self.bam_file = Some(PathBuf::from(content)); },
            Message::SelectBam => { Self::load_file(&mut self.bam_file); },
            Message::MotifChanged(content) => { self.motif_file = Some(PathBuf::from(content)); },
            Message::SelectMotif => { Self::load_file(&mut self.motif_file); },

            Message::OutdirChanged(content) => { self.output = Some(PathBuf::from(content)); },
            Message::SelectOutdir => { Self::load_dir(&mut self.output); },

            Message::RunDante => { self.run(); },
            Message::OpenResults => { self.open_results(); }

            // Message::CheckboxToggled(is_checked) => { self.out_bam = is_checked },
            _ => { todo!(); }
        }
    }

    fn load_file(result: &mut Option<PathBuf>) {
        let path = FileDialog::new().show_open_single_file().unwrap();
        let path = match path {
            Some(path) => path,
            None => return,
        };
        *result = Some(path);
    }

    fn load_dir(result: &mut Option<PathBuf>) {
        let path = FileDialog::new().show_open_single_dir().unwrap();
        let path = match path {
            Some(path) => path,
            None => return,
        };
        *result = Some(path);
    }

    fn run(&self) {
        println!("{:?}", self);

        // required params
        let ref_file: String = self.ref_file.as_ref().unwrap().to_str().unwrap().to_string();
        let bam_file: String = self.bam_file.as_ref().unwrap().to_str().unwrap().to_string();
        let motif_file: String = self.motif_file.as_ref().unwrap().to_str().unwrap().to_string();
        let mut output: String = self.output.as_ref().unwrap().to_str().unwrap().to_string();
        let out_dir = output.clone();
        output.push_str("/remaSTR_result.tsv");

        // optional params
        let out_bam = false;
        let correction = false;
        let dedup = false;
        let print_quality = false;
        let flank = 30;
        let q = 30;
        let score: Option<char> = None;

        run(ref_file, bam_file, motif_file, output.clone(), out_bam, correction, dedup, flank, q, score, print_quality);
        println!("remaSTR finished.");
        // python ../remaSTR_validation/refactoring/src/dante-remaSTR/dante_remastr_standalone.py
        let output_log = Command::new("python")
            .arg("/home/balaz/projects/STRs/remaSTR_validation/refactoring/src/dante-remaSTR/dante_remastr_standalone.py")
            .arg("--verbose")
            .arg("-i")
            .arg(output.clone())
            .arg("--output-dir")
            .arg(out_dir)
            .output().expect("failed to run python part of Dante");
        println!("Dante finished.");
        println!("{:?}", output_log);
    }

    fn open_results(&mut self) {
        let mut output: String = self.output.as_ref().unwrap().to_str().unwrap().to_string();
        output.push_str("/report.html");
        opener::open(output).unwrap();
    }

}

impl State {
    const PAD1: Padding = Padding { left: 0.0, right: 5.0, top: 0.0, bottom: 0.0 };
    const PAD2: Padding = Padding { left: 5.0, right: 0.0, top: 0.0, bottom: 0.0 };
    const LEFT_WIDTH: u16 = 120;

    fn view(&self) -> Element<Message> {
        let tmp: &str = match self.bam_file.as_ref() {
            Some(x) => x.to_str().unwrap(),
            None => "Load some file"
        };
        let pad2 = Padding { left: 5.0, ..Padding::default() };

        column![
            column![
                image("assets/logo_cut.png").height(100),
            ].width(720.0).align_x(Horizontal::Right),
            horizontal_rule(2),
            column![
                Self::loader_row("Reference file:", &self.ref_file, Message::RefChanged, Message::SelectRef),
                Self::loader_row("BAM file:",       &self.bam_file, Message::BamChanged, Message::SelectBam),
                Self::loader_row("Motif file:", &self.motif_file, Message::MotifChanged, Message::SelectMotif),
                horizontal_rule(2),
                Self::loader_row("Output directory:", &self.output, Message::OutdirChanged, Message::SelectOutdir),
                row![
                    container("").width(State::LEFT_WIDTH).padding(State::PAD1),
                    checkbox("Output BAM", self.out_bam).on_toggle(Message::CheckboxToggled),
                    checkbox("Correction", self.correction).on_toggle(Message::CheckboxToggled),
                    checkbox("Dedup", self.dedup).on_toggle(Message::CheckboxToggled),
                    checkbox("Print quality", self.print_quality).on_toggle(Message::CheckboxToggled),
                ].padding(10.0).align_y(Vertical::Center),
                row![
                    container("").width(State::LEFT_WIDTH).padding(State::PAD1),
                    button("Run").on_press(Message::RunDante),
                    text(tmp).size(12),
                ].padding(10.0).align_y(Vertical::Center),
                row![
                    container("").width(State::LEFT_WIDTH).padding(State::PAD1),
                    button("Go to result").style(|_theme: &Theme, _status| {
                        let style = button::Style {
                            background: None,
                            text_color: iced::Color { r: 1.0, g: 1.0, b: 1.0, a: 1.0 },
                            border: iced::Border { color: iced::Color::BLACK, width: 1.0, radius: iced::border::Radius::default() },
                            shadow: iced::Shadow::default()
                        };
                        return style;
                    }).on_press(Message::OpenResults),
                    container(text("Result: ").line_height(2.2).size(12).align_x(Horizontal::Right)),
                ].padding(10.0).align_y(Vertical::Center),
            ].width(720.0).align_x(Horizontal::Left)
        ].into()
    }

    fn loader_row<'a>(desc: &'a str, filename: &'a Option<PathBuf>, on_input: impl Fn(String) -> Message + 'a, on_press: Message) -> Row<'a, Message> {
        let filename_str: &str = match filename.as_ref() {
            Some(x) => x.to_str().unwrap(),
            None => ""
        };

        row![
            container(text(desc).width(State::LEFT_WIDTH).align_x(Horizontal::Right)).padding(State::PAD1),
            text_input("Type file path or click search...", filename_str).on_input(on_input),
            container(button("Search file").on_press(on_press)).padding(State::PAD2),
        ].padding(10.0).align_y(Vertical::Center)
    }
}
