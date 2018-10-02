extern crate run_hulk;

use std::process;

fn main() {
    let config = run_hulk::get_args().expect("Could not get arguments");

    if let Err(e) = run_hulk::run(config) {
        println!("Error: {}", e);
        process::exit(1);
    }
}
