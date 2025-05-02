#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use std::fmt::Display;
use std::fs;
use std::path::Path;

use iced::widget::{column, container, horizontal_rule, image};
use iced::{window, Task};
use iced::window::Settings;
use iced::{Element, Length, Theme, Size, Subscription};

mod welcome_screen;

mod analysis_family;
mod analysis_single;
mod analysis_common;

mod pdf_reporting;
mod metadata_editor;
mod async_tasks;
mod editor_results;

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
    SingleResults(editor_results::Message),
    AnalysisFamily(analysis_family::Message),
    MetadataEditor(metadata_editor::Message),
    Resize(Size),
}

#[derive(Debug)]
enum ContentPage {
    WelcomeScreen(welcome_screen::Data),
    AnalysisSingle(analysis_single::Data),
    SingleResults(editor_results::Data),
    AnalysisFamily(analysis_family::Data),
    MetadataEditor(metadata_editor::Data),
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
            ContentPage::WelcomeScreen(data)  => data.view().map(Message::WelcomeScreen),
            ContentPage::AnalysisSingle(data) => data.view(self.window_size).map(Message::AnalysisSingle),
            ContentPage::AnalysisFamily(data) => data.view(self.window_size).map(Message::AnalysisFamily),
            ContentPage::MetadataEditor(data) => data.view(self.window_size).map(Message::MetadataEditor),
            ContentPage::SingleResults(data)  => data.view(self.window_size).map(Message::SingleResults),
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

    fn update(&mut self, message: Message) -> Task<Message> {
        use ContentPage as CP;
        use metadata_editor::Data as MetaEditor;
        use welcome_screen::analysis_create as analysis_create;
        use welcome_screen::analysis_reopen as analysis_reopen;

        match message {
            // global state changes
            Message::AnalysisSingle(analysis_single::Message::Back) => {
                self.content_page = CP::default(); Task::none()
            }
            Message::AnalysisFamily(analysis_family::Message::Back) => {
                self.content_page = CP::default(); Task::none()
            }
            Message::WelcomeScreen(welcome_screen::Message::CreateAnalysis(name, atype)) => {
                self.content_page = analysis_create(name, atype); Task::none()
            }
            Message::WelcomeScreen(welcome_screen::Message::AnalysisReopen(path)) => {
                self.content_page = analysis_reopen(path); Task::none()
            }

            Message::AnalysisSingle(analysis_single::Message::EditMetadata(source, meta_file)) => {
                let CP::AnalysisSingle(ref data) = self.content_page else { unreachable!() };
                data.save();
                self.content_page = MetaEditor::open(source, meta_file); Task::none()
            }
            Message::AnalysisSingle(analysis_single::Message::EditResults(data)) => {
                data.save();
                self.content_page = editor_results::Data::open(
                    data.get_checked_motif_ids(),
                    data.get_checked_motif_names(),
                    data.get_source(),
                    data.get_sample()
                );
                Task::none()
            }
            Message::MetadataEditor(metadata_editor::Message::Exit(source)) => {
                self.content_page = analysis_reopen(source); Task::none()
            }
            Message::SingleResults(editor_results::Message::Exit(source)) => {
                self.content_page = analysis_reopen(source); Task::none()
            }

            Message::Resize(size) => {
                self.window_size = size; Task::none()
            }

            // relay message
            Message::SingleResults(m) => {
                if let CP::SingleResults(data) = &mut self.content_page {
                    data.update(m);
                };
                Task::none()
            },
            Message::WelcomeScreen(m) => {
                if let CP::WelcomeScreen(data) = &mut self.content_page {
                    data.update(m);
                };
                Task::none()
            },
            Message::AnalysisSingle(m) => {
                if let CP::AnalysisSingle(data) = &mut self.content_page {
                    data.update(m).map(Message::AnalysisSingle)
                } else {
                    Task::none()
                }
            },
            Message::AnalysisFamily(m) => {
                if let CP::AnalysisFamily(data) = &mut self.content_page {
                    data.update(m).map(Message::AnalysisFamily)
                } else {
                    // TODO: kill potential jobs?
                    Task::none()
                }
            },
            Message::MetadataEditor(m) => {
                if let CP::MetadataEditor(data) = &mut self.content_page {
                    data.update(m);
                }
                Task::none()
            },
        }
    }

    fn subscription(&self) -> Subscription<Message> {
        window::resize_events().map(|(_, size)| {
            return Message::Resize(size);
        })
    }
}

fn init_cache<S>(path: &S)
where
    S: AsRef<Path> + Display + ?Sized,
{
    let filenames = [
        "logo.png",
        "STRSet_20220902.tsv",
        "STRSet_20250311.tsv",
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
        include_bytes!("../assets/STRSet_20250311.tsv").to_vec(),
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
    #[allow(non_camel_case_types)]
    STRSet_20250311,
    Custom,
}

impl std::fmt::Display for MotifFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::STRSet_20220902 => "STRSet_20220902",
            Self::STRSet_20250311 => "STRSet_20250311",
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
