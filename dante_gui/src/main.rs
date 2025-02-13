#![windows_subsystem = "windows"]

use iced::alignment::Horizontal;
use iced::font::Weight;
use iced::widget::{column, horizontal_rule, image};
use iced::Font;
use iced::{Element, Padding, Theme};
use native_dialog::FileDialog;
use remastr::run;
use std::ffi::OsStr;
use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;
use std::fmt::Display;

mod welcome_screen;
mod analysis_single;
mod analysis_family;

pub fn main() -> iced::Result {
    let settings = iced::window::Settings {
        size: iced::Size { width: 720.0, height: 480.0 },
        // size: iced::Size{width: 720.0, height: 560.0},
        ..Default::default()
    };
    iced::application("Dante", App::update, App::view)
        .window(settings)
        .theme(|_| Theme::CatppuccinLatte)
        .run()
}

#[derive(Debug, Default)]
struct App {
    content_page: ContentPage,

    bam_file: Option<PathBuf>,
    motif_file: Option<PathBuf>,
    output: Option<PathBuf>,
    out_bam: bool,
    message_line: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ContentPage {
    WelcomeScreen(welcome_screen::Data),
    AnalysisSingle,
    AnalysisFamily,
}

impl Default for ContentPage {
    fn default() -> Self {
        Self::WelcomeScreen(welcome_screen::Data::default())
    }
}

#[derive(Debug, Clone)]
enum Message {
    AnalysisSelected(welcome_screen::Analysis),
    AnalysisCreate,

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

impl App {
    const PAD1: Padding = Padding { left: 0.0, right: 5.0, top: 0.0, bottom: 0.0 };
    const PAD2: Padding = Padding { left: 5.0, right: 0.0, top: 0.0, bottom: 0.0 };
    const LEFT_WIDTH: u16 = 120;
    const BOLD_MONO: Font = Font { weight: Weight::Bold, ..Font::MONOSPACE };
    const DATA_DIR: &str = "dante_data";

    fn view(&self) -> Element<Message> {
        if !Path::new(Self::DATA_DIR).exists() { init_cache(Self::DATA_DIR); }

        let content_area: Element<Message> = match self.content_page {
            ContentPage::WelcomeScreen(data) => welcome_screen::view(self, data),
            ContentPage::AnalysisSingle => analysis_single::view(self),
            ContentPage::AnalysisFamily => analysis_family::view(self),
        };

        column![
            column![
                image(format!("{}/logo.png", Self::DATA_DIR)).height(100)
            ].width(720.0).align_x(Horizontal::Right),
            horizontal_rule(0),
            content_area
        ].into()
    }

    fn update(&mut self, message: Message) {
        match message {
            Message::BamChanged(content) => { self.bam_file = Some(PathBuf::from(content)); },
            Message::SelectBam => { load_file(&mut self.bam_file); },
            Message::MotifChanged(content) => { self.motif_file = Some(PathBuf::from(content)); },
            Message::SelectMotif => { load_file(&mut self.motif_file); },

            Message::OutdirChanged(content) => { self.output = Some(PathBuf::from(content)); },
            Message::SelectOutdir => { load_dir(&mut self.output); },

            Message::RunDante => { run1(self); },
            Message::OpenResults => { open_results(self); },

            Message::CheckboxOutBAM(is_checked) => self.out_bam = is_checked,

            Message::AnalysisSelected(analysis) => { welcome_screen::analysis_set(self, analysis); },
            Message::AnalysisCreate => { welcome_screen::analysis_create(self); }
        }
    }
}

fn load_file(result: &mut Option<PathBuf>) {
    // TODO: "." does not work under Windows
    // let path = FileDialog::new().set_location(".").show_open_single_file().unwrap();
    let path = FileDialog::new().show_open_single_file().unwrap();
    let path = match path {
        Some(path) => path,
        None => return,
    };
    *result = Some(path);
}

fn load_dir(result: &mut Option<PathBuf>) {
    // TODO: "." does not work under Windows
    // let path = FileDialog::new().set_location(".").show_open_single_dir().unwrap();
    let path = FileDialog::new().show_open_single_dir().unwrap();
    let path = match path {
        Some(path) => path,
        None => return,
    };
    *result = Some(path);
}

fn run1(state: &App) {
    println!("{:?}", state);

    // required params
    let Some(ref bam_file) = state.bam_file else {
        return;
    };
    let Some(ref motif_file) = state.motif_file else {
        return;
    };

    let mut output: String = match state.output {
        Some(ref x) => x.display().to_string(),
        None => {
            return;
        },
    };
    let out_dir = output.clone();
    if !Path::new(&out_dir).exists() {
        fs::create_dir(&out_dir).expect("Cannot create directory.");
    }
    output.push_str("/remaSTR_result.tsv");

    // optional params
    let out_bam = state.out_bam;
    let dedup = false;
    let print_quality = false;
    let q = 30;
    let score: Option<char> = None;

    run(bam_file, motif_file, output.clone(), out_bam, (dedup, q, score, print_quality));
    println!("remaSTR finished.");
    // self.message_line = "remaSTR finished.".to_string();
    let bin = format!("{}/dante_remastr_standalone", App::DATA_DIR);
    let output_log = Command::new(bin)
        .arg("--input-tsv").arg(output.clone())
        .arg("--output-dir").arg(out_dir)
        .arg("--verbose")
        .output()
        .expect("failed to run python part of Dante");
    println!("Dante finished.");
    // self.message_line = "Dante finished.".to_string();
    println!("{:?}", output_log);
}

fn open_results(state: &mut App) {
    match state.output.as_ref() {
        Some(x) => {
            let mut output: String = x.to_str().unwrap().to_string();
            output.push_str("/report.html");
            opener::open(output).unwrap();
        },
        None => {
            state.message_line = "No report found.".to_string();
        },
    }
}

fn init_cache<S>(path: &S)
where
    S: AsRef<OsStr> + AsRef<Path> + Display + ?Sized,
{
    let filenames = [
        "logo.png",
        "templates/alignments_template.html",
        "templates/report_template.html",
        "includes/datatables.min.js",
        "includes/jquery-3.6.1.min.js",
        "includes/jquery.dataTables.css",
        "includes/msa.min.gz.js",
        "includes/plotly-2.14.0.min.js",
        "includes/styles.css",
        "includes/w3.css",
    ];

    let contents = [
        include_bytes!("../assets/logo.png").to_vec(),
        include_bytes!("../assets/templates/alignments_template.html").to_vec(),
        include_bytes!("../assets/templates/report_template.html").to_vec(),
        include_bytes!("../assets/includes/datatables.min.js").to_vec(),
        include_bytes!("../assets/includes/jquery-3.6.1.min.js").to_vec(),
        include_bytes!("../assets/includes/jquery.dataTables.css").to_vec(),
        include_bytes!("../assets/includes/msa.min.gz.js").to_vec(),
        include_bytes!("../assets/includes/plotly-2.14.0.min.js").to_vec(),
        include_bytes!("../assets/includes/styles.css").to_vec(),
        include_bytes!("../assets/includes/w3.css").to_vec(),
    ];

    fs::create_dir(path).expect("Cannot create directory.");
    fs::create_dir(format!("{}/templates", path)).expect("Cannot create directory.");
    fs::create_dir(format!("{}/includes", path)).expect("Cannot create directory.");

    for (filename, content) in std::iter::zip(filenames, contents) {
        fs::write(format!("{}/{}", path, filename), content).expect("Unable to write.");
    }

    #[cfg(target_os = "linux")]
    {
        use std::os::unix::fs::PermissionsExt;

        let ctx = include_bytes!("../assets/dante_remastr_standalone").to_vec();
        let bin = format!("{}/dante_remastr_standalone", path);
        fs::write(&bin, ctx).expect("Unable to write.");

        let mut perms = fs::metadata(&bin).unwrap().permissions();
        perms.set_mode(0o700);
        fs::set_permissions(&bin, perms).unwrap();
    }

    #[cfg(target_os = "windows")]
    {
        let ctx = include_bytes!("../assets/dante_remastr_standalone.exe").to_vec();
        let bin = format!("{}/dante_remastr_standalone.exe", path);
        fs::write(bin, ctx).expect("Unable to write.");
    }
}
