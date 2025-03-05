#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use std::fmt::Display;
use std::fs;
use std::path::Path;

use iced::widget::{column, container, horizontal_rule, image};
use iced::window;
use iced::window::Settings;
use iced::{Element, Length, Theme, Size, Subscription};

mod welcome_screen;

mod analysis_family;
mod analysis_single;

mod pdf_reporting;

pub fn main() -> iced::Result {
    if !Path::new(App::DATA_DIR).exists() { init_cache(App::DATA_DIR); }

    iced::application("Dante", App::update, App::view)
        .window(Settings { size: Size { width: 960.0, height: 960.0 }, ..Default::default() })
        .subscription(App::subscription)
        .theme(|_| Theme::CatppuccinLatte)
        .run()
}

#[derive(Debug, Clone)]
enum Message {
    WelcomeScreen(welcome_screen::Message),
    AnalysisSingle(analysis_single::Message),
    AnalysisFamily(analysis_family::Message),
    Resize(Size),
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ContentPage {
    WelcomeScreen(welcome_screen::Data),
    AnalysisSingle(analysis_single::Data),
    AnalysisFamily(analysis_family::Data),
}

impl Default for ContentPage {
    fn default() -> Self {
        Self::WelcomeScreen(welcome_screen::Data::default())
    }
}

#[derive(Debug, Default)]
struct App {
    window_size: Size,
    content_page: ContentPage,
}

use iced::Padding;
use serde::{Deserialize, Serialize};
impl App {
    const PAD1: Padding = Padding { left: 0.0, right: 5.0, top: 0.0, bottom: 0.0 };
    const PAD2: Padding = Padding { left: 5.0, right: 0.0, top: 0.0, bottom: 0.0 };
    const DATA_DIR: &str = "dante_data";
    const H1_SIZE: u16 = 26;

    fn view(&self) -> Element<Message> {
        let content_area: Element<Message> = match &self.content_page {
            ContentPage::WelcomeScreen(data) => data.view().map(Message::WelcomeScreen),
            ContentPage::AnalysisSingle(data) => data.view(self.window_size).map(Message::AnalysisSingle),
            ContentPage::AnalysisFamily(data) => data.view(self.window_size).map(Message::AnalysisFamily),
        };

        // let content_area = std::convert::Into::<Element<Message>>::into(content_area).explain(iced::Color::BLACK);
        use iced::alignment::Horizontal;
        use iced::widget::container::background;
        use iced::Color;
        column![
            container(image(format!("{}/logo.png", Self::DATA_DIR)).width(900).height(125))
                .width(Length::Fill).align_x(Horizontal::Center)
                .style(|_| { background(Color { r: 0.77, g: 0.82, b: 0.84, a: 1.0 }) }),
            horizontal_rule(1),
            content_area
        ].align_x(Horizontal::Center).into()
    }

    fn update(&mut self, message: Message) {
        use ContentPage as CP;
        match message {
            // global state changes
            Message::WelcomeScreen(welcome_screen::Message::CreateAnalysis)
                => { welcome_screen::analysis_create(self); }
            Message::WelcomeScreen(welcome_screen::Message::AnalysisReopen(path))
                => { welcome_screen::analysis_reopen(self, path); }
            Message::AnalysisSingle(analysis_single::Message::Back)
                => { back(self); }
            Message::AnalysisFamily(analysis_family::Message::Back)
                => { back(self); }
            Message::Resize(size)
                => { self.window_size = size; }

            // local state changes
            Message::WelcomeScreen(m) => {
                if let CP::WelcomeScreen(data) = &mut self.content_page {
                    data.update(m);
                }
            },
            Message::AnalysisSingle(m) => {
                if let CP::AnalysisSingle(data) = &mut self.content_page {
                    data.update(m);
                }
            },
            Message::AnalysisFamily(m) => {
                if let CP::AnalysisFamily(data) = &mut self.content_page {
                    data.update(m);
                }
            },
        }
    }

    fn subscription(&self) -> Subscription<Message> {
        window::resize_events().map(|(_, size)| {
            return Message::Resize(size);
        })
    }
}

fn back(state: &mut App) {
    state.content_page = ContentPage::WelcomeScreen(welcome_screen::Data::default());
}

fn init_cache<S>(path: &S)
where
    S: AsRef<Path> + Display + ?Sized,
{
    let filenames = [
        "logo.png",
        "STRSet_20220902.tsv",
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
        include_bytes!("../assets/STRSet_20220902.tsv").to_vec(),
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
enum MotifFile {
    #[allow(non_camel_case_types)]
    STRSet_20220902,
    Custom,
}

impl std::fmt::Display for MotifFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::STRSet_20220902 => "STRSet_20220902",
            Self::Custom => "custom",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
enum Sex {
    Male,
    Female,
    Unknown,
}

impl std::fmt::Display for Sex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Male => "male",
            Self::Female => "female",
            Self::Unknown => "unknown",
        })
    }
}
