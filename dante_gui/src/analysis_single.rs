use iced::alignment::{Horizontal, Vertical};
use iced::widget::{button, checkbox, column, container, horizontal_rule, row, text, text_input};
use iced::Element;
use std::path::PathBuf;
use std::env;
use std::path::Path;

use crate::App;
use crate::Message;

pub fn view(state: &App) -> Element<Message> {
    column![
        button("Back").on_press(Message::Back),
        loader_row("BAM file:",       &state.bam_file, Message::BamChanged, Message::SelectBam),
        loader_row("Motif file:", &state.motif_file, Message::MotifChanged, Message::SelectMotif),
        horizontal_rule(2),
        loader_row("Output directory:", &state.output, Message::OutdirChanged, Message::SelectOutdir),

        row![
            container("").width(App::LEFT_WIDTH).padding(App::PAD1),
            checkbox("Output BAM", state.out_bam).on_toggle(Message::CheckboxOutBAM),
        ].padding(10.0).align_y(Vertical::Center),

        run_button(state),
        draw_open_button(state),
    ].width(720.0).align_x(Horizontal::Left).into()
}

fn run_button<'a>(state: &App) -> Element<'a, Message> {
    row![
        container("").width(App::LEFT_WIDTH).padding(App::PAD1),
        button("Run").on_press(Message::RunDante),
        container(text(state.message_line.clone()).align_x(Horizontal::Left)).padding(App::PAD2),
    ].padding(10.0).align_y(Vertical::Center).into()
}

fn draw_open_button<'a>(state: &App) -> Element<'a, Message> {
    let report_present;
    let report_line;
    match &state.output {
        Some(x) => {
            let mut x = path_to_string(x);
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
            container("").width(App::LEFT_WIDTH).padding(App::PAD1),
            button("Open results").on_press(Message::OpenResults),
            container(text(report_line).align_x(Horizontal::Left)).padding(App::PAD2),
        ].padding(10.0).align_y(Vertical::Center).into()
    } else {
        row![
            container("").width(App::LEFT_WIDTH).padding(App::PAD1),
            button("Open results"),
            container(text(report_line).align_x(Horizontal::Left)).padding(App::PAD2),
        ].padding(10.0).align_y(Vertical::Center).into()
    }
}

fn loader_row<'a>(
    desc: &'a str, filename: &'a Option<PathBuf>, on_input: impl Fn(String) -> Message + 'a, on_press: Message
) -> Element<'a, Message> {
    let filename_str: String = match filename.as_ref() {
        Some(x) => path_to_string(x),
        None => "".to_string()
    };

    row![
        container(text(desc).width(App::LEFT_WIDTH).align_x(Horizontal::Right)).padding(App::PAD1),
        text_input("Type path or click search...", &filename_str).on_input(on_input).font(App::BOLD_MONO),
        container(button("Search").on_press(on_press)).padding(App::PAD2),
    ].padding(10.0).align_y(Vertical::Center).into()
}

fn path_to_string(path: &Path) -> String {
    let cwd = env::current_dir().unwrap().display().to_string();
    match path.strip_prefix(cwd) {
        Ok(x) => { x.display().to_string() },
        Err(_) => { path.display().to_string() }
    }
}


