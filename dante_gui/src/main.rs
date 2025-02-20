#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use iced::Element;
use iced::widget::{column, horizontal_rule, image};
use std::fs;
use std::path::Path;
use std::fmt::Display;

mod welcome_screen;
mod analysis_single;
mod analysis_family;

pub fn main() -> iced::Result {
    if !Path::new(App::DATA_DIR).exists() { init_cache(App::DATA_DIR); }
    iced::application("Dante", App::update, App::view).theme(|_| iced::Theme::CatppuccinLatte).run()
}

#[derive(Debug, Clone)]
enum Message {
    WelcomeScreen(welcome_screen::Message),
    AnalysisSingle(analysis_single::Message),
    AnalysisFamily(analysis_family::Message),
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
    content_page: ContentPage,
}

use iced::{Padding, Font};
impl App {
    const PAD1: Padding = Padding { left: 0.0, right: 5.0, top: 0.0, bottom: 0.0 };
    const PAD2: Padding = Padding { left: 5.0, right: 0.0, top: 0.0, bottom: 0.0 };
    const LEFT_WIDTH: u16 = 120;
    const BOLD_MONO: Font = Font { weight: iced::font::Weight::Bold, ..Font::MONOSPACE };
    const DATA_DIR: &str = "dante_data";

    fn view(&self) -> Element<Message> {
        let content_area: Element<Message> = match &self.content_page {
            ContentPage::WelcomeScreen(data) => welcome_screen::view(data).map(Message::WelcomeScreen),
            ContentPage::AnalysisSingle(data) => analysis_single::view(data).map(Message::AnalysisSingle),
            ContentPage::AnalysisFamily(data) => data.view().map(Message::AnalysisFamily),
        };

        // let content_area = std::convert::Into::<Element<Message>>::into(content_area).explain(iced::Color::BLACK);
        column![
            image(format!("{}/logo.png", Self::DATA_DIR)).width(900).height(125),
            horizontal_rule(0),
            content_area
        ].align_x(iced::alignment::Horizontal::Center).into()
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

            // local state changes
            Message::WelcomeScreen(m) => {
                if let CP::WelcomeScreen(data) = &mut self.content_page {
                    welcome_screen::update(data, m); 
                }
            }
            Message::AnalysisSingle(m) => {
                if let CP::AnalysisSingle(data) = &mut self.content_page {
                    analysis_single::update(data, m);
                }
            }
            Message::AnalysisFamily(m) => {
                if let CP::AnalysisFamily(data) = &mut self.content_page {
                    data.update(m);
                }
            }
        }
    }
}

fn back(state: &mut App){
    state.content_page = ContentPage::WelcomeScreen(welcome_screen::Data::default());
}

fn init_cache<S>(path: &S)
where
    S: AsRef<Path> + Display + ?Sized,
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
