#![cfg_attr(all(windows, not(debug_assertions)), windows_subsystem = "windows")]

fn main() {
    malaria_amplicon_nf_desktop_lib::run()
}
