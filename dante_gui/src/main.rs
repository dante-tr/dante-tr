#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use std::fs;
use std::path::{Path, PathBuf};

use iced::widget::{column, container, horizontal_rule, image, row, text};
use iced::{window, Task};
use iced::window::Settings;
use iced::{Element, Length, Theme, Size, Subscription};
use iced::Padding;
use serde::{Deserialize, Serialize};

mod welcome_screen;

mod analysis_family;
mod analysis_single;
mod analysis_common;

mod pdf_reporting;
mod metadata_editor;
mod async_tasks;
mod editor_results;

mod components;

// defines const EMBEDDED_FILES: [(&str, &[u8]); N];
include!(concat!(env!("OUT_DIR"), "/embedded_assets.rs"));

pub fn main() -> iced::Result {
    if !Path::new(App::DATA_DIR).exists() { App::init_cache(); }

    iced::application("Dante", App::update, App::view)
        .window(Settings { size: Size { width: 960.0, height: 960.0 }, ..Default::default() })
        .subscription(App::subscription)
        .theme(|_| Theme::CatppuccinLatte)
        .run()
}

#[derive(Debug, Clone)]
#[allow(clippy::large_enum_variant)]
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

impl App {
    const PAD1: Padding = Padding { left: 0.0, right: 5.0, top: 0.0, bottom: 0.0 };
    const PAD2: Padding = Padding { left: 5.0, right: 0.0, top: 0.0, bottom: 0.0 };
    const DATA_DIR: &str = "dante_data";
    const H1_SIZE: u16 = 26;

    fn view(&self) -> Element<'_, Message> {
        let header: Element<Message> = self.view_header();
        let content_area: Element<Message> = match &self.content_page {
            ContentPage::WelcomeScreen(data)  => data.view().map(Message::WelcomeScreen),
            ContentPage::AnalysisSingle(data) => data.view(self.window_size).map(Message::AnalysisSingle),
            ContentPage::AnalysisFamily(data) => data.view(self.window_size).map(Message::AnalysisFamily),
            ContentPage::MetadataEditor(data) => data.view(self.window_size).map(Message::MetadataEditor),
            ContentPage::SingleResults(data)  => data.view(self.window_size).map(Message::SingleResults),
        };

        // let content_area = std::convert::Into::<Element<Message>>::into(content_area).explain(iced::Color::BLACK);
        column![
            header,
            horizontal_rule(1),
            content_area,
        ].into()
    }

    fn view_header(&self) -> Element<'_, Message> {
        use iced::alignment::Horizontal;
        use iced::widget::container::background;
        use iced::Color;
        let version = format!("v{} ", env!("CARGO_PKG_VERSION"));
        let logo = Self::get_filename("assets/includes/logo.png");
        column![
            row![
                container(text(version)).width(Length::Fill).align_x(Horizontal::Right)
                    .style(|_| { background(Color { r: 0.77, g: 0.82, b: 0.84, a: 1.0 }) })
            ],
            container(image(logo).width(900).height(125))
                .width(Length::Fill).align_x(Horizontal::Center)
                .style(|_| { background(Color { r: 0.77, g: 0.82, b: 0.84, a: 1.0 }) }),
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
                let CP::AnalysisSingle(ref mut data) = self.content_page else { unreachable!() };
                data.save();
                self.content_page = MetaEditor::open(source, meta_file); Task::none()
            }
            Message::AnalysisSingle(analysis_single::Message::EditResults(mut data)) => {
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

    fn get_filename(old_filename: &str) -> PathBuf {
        let tmp = PathBuf::from(old_filename);
        let new_filepath = PathBuf::from(
            format!("{}/{}", Self::DATA_DIR, tmp.strip_prefix("assets").unwrap().display())
        );
        return new_filepath;
    }

    fn init_cache() {
        // set correct destinations
        let mut table: Vec<(PathBuf, &[u8])> = Vec::new();
        for (old_filename, content) in EMBEDDED_FILES {
            let new_filepath = Self::get_filename(old_filename);
            // println!("{} -> {}", old_filename, new_filepath.display());
            table.push((new_filepath, content));
        }

        // create directories
        for (filepath, _) in &table {
            let tmp = filepath.parent().unwrap();
            match fs::create_dir_all(tmp) {  // This works, but strictly following the documentation shouldn't
                Ok(_)  => { /* println!("Creating {}: Ok.", tmp.display()) */ },
                Err(_) => { println!("Creating {}: Failed!", tmp.display()) }
            }
        }

        for (filepath, content) in &table {
            fs::write(filepath, content).expect("Unable to write.");
        }

        println!("Assets extracted to {}/", Self::DATA_DIR);
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
