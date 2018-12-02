extern crate clap;
extern crate csv;
extern crate regex;
extern crate walkdir;

use clap::{App, Arg};
use regex::Regex;
use std::collections::HashMap;
use std::error::Error;
use std::process::{Command, Stdio};
use std::{
    env,
    ffi::OsStr,
    fs::{self, DirBuilder, File},
    io::{Read, Write},
    path::{Path, PathBuf},
};
use walkdir::WalkDir;

// --------------------------------------------------
// Custom types
// --------------------------------------------------
type Record = HashMap<String, String>;

type MyResult<T> = Result<T, Box<Error>>;

#[derive(Debug)]
pub struct Config {
    alias_file: Option<String>,
    bin_dir: Option<String>,
    kmer_size: Option<u32>,
    min_kmer_count: Option<u32>,
    interval: Option<u32>,
    sketch_size: Option<u32>,
    num_threads: Option<u32>,
    out_dir: PathBuf,
    query: Vec<String>,
    reads_are_fasta: bool,
    create_weighted_matrix: bool,
}

// --------------------------------------------------
pub fn run(config: Config) -> MyResult<()> {
    println!("config {:?}", config);

    let files = find_files(&config.query, None)?;

    if files.len() == 0 {
        let msg = format!("No input files from query \"{:?}\"", &config.query);
        return Err(From::from(msg));
    }

    println!(
        "Will process {} file{}",
        files.len(),
        if files.len() == 1 { "" } else { "s" }
    );

    let out_dir = &config.out_dir;
    if !out_dir.is_dir() {
        DirBuilder::new().recursive(true).create(&out_dir)?;
    }

    let sketch_dir = sketch_files(&config, &files)?;
    println!("Sketch dir = {:?}", sketch_dir);

    let fig_dir = smash_sketches(&config, &sketch_dir)?;
    println!("Done, see output in \"{}\"", fig_dir.to_string_lossy());

    Ok(())
}

// --------------------------------------------------
pub fn get_args() -> MyResult<Config> {
    let matches = App::new("HULK")
        .version("0.1.0")
        .author("Ken Youens-Clark <kyclark@email.arizona.edu>")
        .about("Run HULK")
        .arg(
            Arg::with_name("query")
                .short("q")
                .long("query")
                .value_name("FILE_OR_DIR")
                .help("File or input directory")
                .required(true)
                .min_values(1),
        )
        .arg(
            Arg::with_name("out_dir")
                .short("o")
                .long("out_dir")
                .value_name("DIR")
                .help("Output directory"),
        )
        .arg(
            Arg::with_name("alias")
                .short("a")
                .long("alias")
                .value_name("FILE")
                .help("Aliases for sample names"),
        )
        .arg(
            Arg::with_name("kmer_size")
                .short("k")
                .long("kmer_size")
                .value_name("INT")
                .default_value("11")
                .help("K-mer size"),
        )
        .arg(
            Arg::with_name("min_kmer_count")
                .short("m")
                .long("min_kmer_count")
                .value_name("INT")
                .default_value("1")
                .help("Minimum k-mer count"),
        )
        .arg(
            Arg::with_name("interval")
                .short("i")
                .long("interval")
                .value_name("INT")
                .default_value("0")
                .help("Size of read sampling interval"),
        )
        .arg(
            Arg::with_name("sketch_size")
                .short("s")
                .long("sketch_size")
                .value_name("INT")
                .default_value("256")
                .help("Sketch size"),
        )
        .arg(
            Arg::with_name("num_threads")
                .short("t")
                .long("num_threads")
                .value_name("INT")
                .default_value("8")
                .help("Number of threads"),
        )
        .arg(
            Arg::with_name("reads_are_fasta")
                .short("f")
                .long("reads_are_fasta")
                .help("Input reads are in FASTA format"),
        )
        .arg(
            Arg::with_name("create_weighted_matrix")
                .short("w")
                .long("create_weighted_matrix")
                .help("Create a pairwise weighted Jaccard Similarity matrix"),
        )
        .arg(
            Arg::with_name("bin_dir")
                .short("b")
                .long("bin_dir")
                .value_name("DIR")
                .help("Location of binaries"),
        )
        .get_matches();

    let out_dir = match matches.value_of("out_dir") {
        Some(x) => PathBuf::from(x),
        _ => {
            let cwd = env::current_dir()?;
            cwd.join(PathBuf::from("hulk-out"))
        }
    };

    let alias = matches.value_of("alias").and_then(|x| Some(x.to_string()));

    let num_threads = matches
        .value_of("num_threads")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let kmer_size = matches
        .value_of("kmer_size")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let min_kmer_count = matches
        .value_of("min_kmer_count")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let interval = matches
        .value_of("min_kmer_count")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let bin_dir = match matches.value_of("bin_dir") {
        Some(x) => Some(x.to_string()),
        _ => None,
    };

    let sketch_size = matches
        .value_of("sketch_size")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let config = Config {
        alias_file: alias,
        bin_dir: bin_dir,
        kmer_size: kmer_size,
        min_kmer_count: min_kmer_count,
        interval: interval,
        sketch_size: sketch_size,
        num_threads: num_threads,
        out_dir: out_dir,
        query: matches.values_of_lossy("query").unwrap(),
        reads_are_fasta: matches.is_present("reads_are_fasta"),
        create_weighted_matrix: matches.is_present("create_weighted_matrix"),
    };

    Ok(config)
}

// --------------------------------------------------
fn find_files(paths: &Vec<String>, re: Option<&Regex>) -> MyResult<Vec<String>> {
    let mut files = vec![];
    for path in paths {
        let meta = fs::metadata(path)?;
        if meta.is_file() {
            files.push(path.to_owned());
        } else {
            for entry in fs::read_dir(path)? {
                let entry = entry?;
                let meta = entry.metadata()?;
                let name = entry.path().display().to_string();
                if meta.is_file() {
                    let ok = match re {
                        Some(r) => r.is_match(&name),
                        _ => true,
                    };

                    if ok {
                        files.push(name);
                    }
                }
            }
        };
    }

    files.sort();
    Ok(files)
}

// --------------------------------------------------
fn sketch_files(config: &Config, files: &Vec<String>) -> MyResult<PathBuf> {
    let sketch_dir = config.out_dir.join(PathBuf::from("sketches"));
    if !sketch_dir.is_dir() {
        DirBuilder::new().recursive(true).create(&sketch_dir)?;
    }

    let mut args = vec![];

    if let Some(t) = config.num_threads {
        if t > 0 && t < 64 {
            args.push(format!("-p {}", t));
        }
    }

    if let Some(k) = config.kmer_size {
        args.push(format!("-k {}", k));
    }

    if let Some(s) = config.sketch_size {
        args.push(format!("-s {}", s));
    }

    if let Some(i) = config.interval {
        if i > 0 {
            args.push(format!("-i {}", i));
        }
    }

    if config.reads_are_fasta {
        args.push("--fasta".to_string());
    }

    let aliases = get_aliases(&config.alias_file)?;
    let mut jobs = vec![];

    for file in files.iter() {
        let basename = basename(&file, &aliases);
        let out_file = sketch_dir.join(basename);
        let hulk_file = format!("{}.sketch", out_file.display());

        if !Path::new(&hulk_file).exists() {
            jobs.push(format!(
                "hulk sketch {} -o {} -f {}",
                args.join(" "),
                out_file.display(),
                file,
            ));
        }
    }

    if jobs.len() > 0 {
        run_jobs(&jobs, "Sketching files", 8)?;
    } else {
        println!("No sketch jobs to run, skipping this step");
    }

    let re = Regex::new(r"\.sketch$").unwrap();
    let mut sketches: Vec<String> = WalkDir::new(&sketch_dir)
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| !e.file_type().is_dir())
        .map(|e| e.path().display().to_string())
        .filter(|name| re.is_match(name))
        .collect();

    sketches.sort();

    if files.len() != sketches.len() {
        return Err(From::from("Failed to create all sketches"));
    }

    Ok(sketch_dir)
}

// --------------------------------------------------
fn run_jobs(jobs: &Vec<String>, msg: &str, num_concurrent: u32) -> MyResult<()> {
    let num_jobs = jobs.len();

    if num_jobs > 0 {
        println!(
            "{} (# {} job{} @ {})",
            msg,
            num_jobs,
            if num_jobs == 1 { "" } else { "s" },
            num_concurrent
        );

        let mut process = Command::new("parallel")
            .arg("-j")
            .arg(num_concurrent.to_string())
            .arg("--halt")
            .arg("soon,fail=1")
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .spawn()?;

        {
            let stdin = process.stdin.as_mut().expect("Failed to open stdin");
            stdin
                .write_all(jobs.join("\n").as_bytes())
                .expect("Failed to write to stdin");
        }

        let result = process.wait()?;
        if !result.success() {
            return Err(From::from("Failed to run jobs in parallel"));
        }
    }

    Ok(())
}

// --------------------------------------------------
fn get_aliases(alias_file: &Option<String>) -> Result<Option<Record>, Box<Error>> {
    match alias_file {
        None => Ok(None),
        Some(file) => {
            let alias_fh = match File::open(file) {
                Ok(file) => file,
                Err(e) => {
                    let msg = format!("Failed to open \"{}\": {}", file, e.to_string());
                    return Err(From::from(msg));
                }
            };

            let mut aliases = HashMap::new();
            let delimiter = match Path::new(&file).extension() {
                Some(ext) => match ext.to_str() {
                    Some("csv") => b',',
                    _ => b'\t',
                },
                _ => b'\t',
            };

            let mut rdr = csv::ReaderBuilder::new()
                .delimiter(delimiter)
                .from_reader(alias_fh);

            for result in rdr.deserialize() {
                let record: Record = result?;
                let name = record.get("sample_name");
                let alias = record.get("alias");

                match (name, alias) {
                    (Some(name), Some(alias)) => {
                        aliases.insert(name.to_string(), alias.to_string());
                        ()
                    }
                    _ => println!("Missing sample_name or alias"),
                }
            }

            if aliases.len() > 0 {
                Ok(Some(aliases))
            } else {
                Ok(None)
            }
        }
    }
}

// --------------------------------------------------
fn smash_sketches(config: &Config, sketch_dir: &PathBuf) -> MyResult<PathBuf> {
    let matrix_arg = format!(
        "--{}jsMatrix",
        if config.create_weighted_matrix {
            "w"
        } else {
            ""
        }
    );

    let fig_dir = config.out_dir.join("figures");
    if !fig_dir.is_dir() {
        DirBuilder::new().recursive(true).create(&fig_dir)?;
    }

    let out_prefix = "hulk";
    let smash_out = fig_dir.join(&out_prefix);

    let hulk_out = Command::new("hulk")
        .arg("smash")
        .arg(matrix_arg)
        .arg("-d")
        .arg(format!("{}/", sketch_dir.to_string_lossy()))
        .arg("--outFile")
        .arg(&smash_out)
        .output()?;

    if !hulk_out.status.success() {
        return Err(From::from(format!(
            "Error ({}): {}",
            hulk_out.status,
            String::from_utf8_lossy(&hulk_out.stderr)
        )));
    }

    let out_file_name = format!(
        "{}.{}js-matrix.csv",
        out_prefix,
        if config.create_weighted_matrix {
            "w"
        } else {
            ""
        }
    );

    let out_file = fig_dir.join(out_file_name);
    if !out_file.is_file() {
        let msg = format!(
            "Failed to create HULK output \"{}\"",
            out_file.to_string_lossy()
        );
        return Err(From::from(msg));
    }

    let mut f = File::open(&out_file)?;
    let mut contents = String::new();
    f.read_to_string(&mut contents)?;
    let lines: Vec<&str> = contents.split("\n").collect();
    let file_names: Vec<&str> = lines[0]
        .split(",")
        .map(|x| Path::new(basename(&x, &None)))
        .map(|p| p.file_stem().and_then(OsStr::to_str).unwrap())
        .collect();

    let dist = fig_dir.join("distance.tab");
    let mut out = File::create(&dist)?;

    for (i, line) in lines.iter().enumerate() {
        if line == &"" {
            break;
        } else if i == 0 {
            write!(out, "\t{}\n", file_names.join("\t"))?;
        } else {
            let n: Vec<String> = line
                .split(",")
                .map(|s| s.parse::<f64>())
                .filter_map(Result::ok)
                .map(|n| format!("{:.04}", (1. - (n / 100.0)).to_string()))
                .collect();

            write!(out, "{}\t{}\n", &file_names[i - 1], n.join("\t"))?;
        }
    }

    let make_figures = "make_figures.r";
    let make_figures_path = match &config.bin_dir {
        Some(d) => Path::new(&d).join(make_figures),
        _ => PathBuf::from(make_figures),
    };

    println!("Making figures");
    let mk_figs = Command::new(&make_figures_path)
        .arg("-o")
        .arg(&fig_dir)
        .arg("-m")
        .arg(&dist)
        .output();

    if let Err(e) = mk_figs {
        let msg = format!(
            "Failed to run \"{}\": {}",
            make_figures_path.to_string_lossy(),
            e.to_string()
        );
        return Err(From::from(msg));
    };

    Ok(fig_dir)
}

// --------------------------------------------------
fn basename<'a>(filename: &'a str, aliases: &'a Option<Record>) -> &'a str {
    let mut parts: Vec<&str> = filename.split("/").collect();
    let name = match parts.pop() {
        Some(x) => x,
        None => filename,
    };

    if let Some(a) = aliases {
        match a.get(name) {
            Some(alias) => alias,
            _ => name,
        }
    } else {
        name
    }
}
