//! This program identifies potential repeat nodes in a GFA graph.

use anyhow::Result;
use clap::{arg, command, value_parser, ArgMatches};
use std::path::PathBuf;

mod gaf;
mod gfa;

fn cli() -> ArgMatches {
    command!()
        .about("A Bidirected Repeat Path Enumerator")
        .arg(arg!(<GFA> "Input file in GFA format.").value_parser(value_parser!(PathBuf)))
        .arg(
            arg!(-r --repeat [REPEAT] "Repeat node size limit")
                .value_parser(value_parser!(usize))
                .default_value("10000"),
        )
        .arg(
            arg!(-n --neighbor [NEIGHBOR] "Minimum neighboring node size")
                .value_parser(value_parser!(usize))
                .default_value("10000"),
        )
        .arg(
            arg!(-i --inout [INOUT] "In/out degree threshold")
                .value_parser(value_parser!(usize))
                .default_value("2"),
        )
        .arg(
            arg!(-g --gaf <GAF> "Input GAF file from GraphAligner.")
                .value_parser(value_parser!(PathBuf)),
        )
        .get_matches()
}

fn main() -> Result<()> {
    let args = cli();
    let gaf = args.get_one::<PathBuf>("gaf").cloned();

    // print nodes
    let nodes = gfa::nodes(args, gaf.is_some())?;

    // optionally print paths from the GAF
    match nodes {
        Some(n) => gaf::count_gaf_paths(gaf.unwrap(), n),
        // end here
        None => Ok(()),
    }
}
