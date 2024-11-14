use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, checkbox, column, container, horizontal_rule, image, row, text, text_input, Row};
use iced::{Element, Padding, Theme};
use native_dialog::FileDialog;
use std::path::PathBuf;
use remastr::run;
use std::process::Command;
use std::env;
use std::path::Path;
use iced::Font;
use iced::font::Weight;

pub fn main() -> iced::Result {
    let settings = iced::window::Settings {
        size: iced::Size{width: 720.0, height: 480.0},
        ..Default::default()
    };
    iced::application("Dante", State::update, State::view)
        .window(settings)
        .theme(|_| { Theme::CatppuccinLatte })
        .run()
}

#[derive(Debug, Default)]
struct State {
    // required params
    ref_file: Option<PathBuf>,
    bam_file: Option<PathBuf>,
    motif_file: Option<PathBuf>,
    output: Option<PathBuf>,

    // optional params
    out_bam: bool,
    // correction: bool, dedup: bool, _flank: usize, _q: u8, _score: Option<char>, print_quality: bool,

    // GUI specifics
    message_line: String,
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
    CheckboxOutBAM(bool),
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

            Message::CheckboxOutBAM(is_checked) => { self.out_bam = is_checked },
        }
    }

    fn load_file(result: &mut Option<PathBuf>) {
        let path = FileDialog::new().set_location(".").show_open_single_file().unwrap();
        let path = match path {
            Some(path) => path,
            None => return,
        };
        *result = Some(path);
    }

    fn load_dir(result: &mut Option<PathBuf>) {
        let path = FileDialog::new().set_location(".").show_open_single_dir().unwrap();
        let path = match path {
            Some(path) => path,
            None => return,
        };
        *result = Some(path);
    }

    fn run(&self) {
        println!("{:?}", self);

        // required params
        let Some(ref ref_file) = self.ref_file else { return; };
        let Some(ref bam_file) = self.bam_file else { return; };
        let Some(ref motif_file) = self.motif_file else { return; };

        let mut output: String = match self.output {
            Some(ref x) => { x.display().to_string() },
            None => { return; }
        };
        let out_dir = output.clone();
        output.push_str("/remaSTR_result.tsv");

        // optional params
        let out_bam = self.out_bam;
        let correction = false;
        let dedup = false;
        let print_quality = false;
        let flank = 30;
        let q = 30;
        let score: Option<char> = None;

        run(ref_file, bam_file, motif_file, output.clone(), out_bam, correction, dedup, flank, q, score, print_quality);
        println!("remaSTR finished.");
        // self.message_line = "remaSTR finished.".to_string();
        let output_log = Command::new("python")
            .arg("/home/balaz/projects/STRs/remaSTR_validation/refactoring/src/dante-remaSTR/dante_remastr_standalone.py")
            .arg("--input-tsv").arg(output.clone())
            .arg("--output-dir").arg(out_dir)
            .arg("--verbose")
            .output().expect("failed to run python part of Dante");
        println!("Dante finished.");
        // self.message_line = "Dante finished.".to_string();
        println!("{:?}", output_log);
    }

    fn open_results(&mut self) {
        match self.output.as_ref() {
            Some(x) => { 
                let mut output: String = x.to_str().unwrap().to_string();
                output.push_str("/report.html");
                opener::open(output).unwrap();
            }
            None => {
                self.message_line = "No report found.".to_string();
            }
        }
    }

}

impl State {
    const PAD1: Padding = Padding { left: 0.0, right: 5.0, top: 0.0, bottom: 0.0 };
    const PAD2: Padding = Padding { left: 5.0, right: 0.0, top: 0.0, bottom: 0.0 };
    const LEFT_WIDTH: u16 = 120;
    const BOLD_MONO: Font = Font { weight: Weight::Bold, ..Font::MONOSPACE };

    fn view(&self) -> Element<Message> {
        column![
            column![
                image("assets/logo_v3.png").height(100),
            ].width(720.0).align_x(Horizontal::Right),
            horizontal_rule(0),
            column![
                Self::loader_row("Reference file:", &self.ref_file, Message::RefChanged, Message::SelectRef),
                Self::loader_row("BAM file:",       &self.bam_file, Message::BamChanged, Message::SelectBam),
                Self::loader_row("Motif file:", &self.motif_file, Message::MotifChanged, Message::SelectMotif),
                horizontal_rule(2),
                Self::loader_row("Output directory:", &self.output, Message::OutdirChanged, Message::SelectOutdir),

                row![
                    container("").width(State::LEFT_WIDTH).padding(State::PAD1),
                    checkbox("Output BAM", self.out_bam).on_toggle(Message::CheckboxOutBAM),
                ].padding(10.0).align_y(Vertical::Center),

                self.run_button(),
                self.draw_open_button(),
            ].width(720.0).align_x(Horizontal::Left)
        ].into()
    }

    fn loader_row<'a>(desc: &'a str, filename: &'a Option<PathBuf>, on_input: impl Fn(String) -> Message + 'a, on_press: Message) -> Row<'a, Message> {
        let filename_str: String = match filename.as_ref() {
            Some(x) => Self::path_to_string(x),
            None => "".to_string()
        };

        row![
            container(text(desc).width(State::LEFT_WIDTH).align_x(Horizontal::Right)).padding(State::PAD1),
            text_input("Type path or click search...", &filename_str).on_input(on_input).font(State::BOLD_MONO),
            container(button("Search").on_press(on_press)).padding(State::PAD2),
        ].padding(10.0).align_y(Vertical::Center)
    }

    fn run_button<'a>(&self) -> Row<'a, Message> {
        row![
            container("").width(State::LEFT_WIDTH).padding(State::PAD1),
            button("Run").on_press(Message::RunDante),
            container(text(self.message_line.clone()).align_x(Horizontal::Left)).padding(State::PAD2),
        ].padding(10.0).align_y(Vertical::Center)
    }

    fn draw_open_button<'a>(&self) -> Row<'a, Message> {
        let report_present;
        let report_line;
        match &self.output {
            Some(x) => {
                let mut x = Self::path_to_string(x);
                x.push_str("/report.html");
                if Path::new(&x).exists() {
                    report_present = true;
                    report_line = format!("Report file stored in {}.", x);
                } else {
                    report_present = false;
                    report_line = "No report file present.".to_string();
                }
                // Path::new("/etc/hosts").exists()
            },
            None => {
                report_present = false;
                report_line = "No report file present.".to_string();
            }
        };

        if report_present {
            row![
                container("").width(State::LEFT_WIDTH).padding(State::PAD1),
                button("Open results").on_press(Message::OpenResults),
                container(text(report_line).align_x(Horizontal::Left)).padding(State::PAD2),
            ].padding(10.0).align_y(Vertical::Center)
        } else {
            row![
                container("").width(State::LEFT_WIDTH).padding(State::PAD1),
                button("Open results"),
                container(text(report_line).align_x(Horizontal::Left)).padding(State::PAD2),
            ].padding(10.0).align_y(Vertical::Center)
        }
    }

    fn path_to_string(path: &Path) -> String {
        let cwd = env::current_dir().unwrap().display().to_string();
        match path.strip_prefix(cwd) {
            Ok(x) => { x.display().to_string() },
            Err(_) => { path.display().to_string() }
        }
    }
}
