//! Bidirected Path Enumerator with Flexible Orientation Traversal
//! This program identifies potential repeat nodes in a GFA graph.

use anyhow::{Context, Result};
use clap::{arg, command, value_parser, ArgMatches};
use gfa::gfa::{Orientation, GFA};
use gfa::parser::GFAParser;
use std::collections::HashMap;
use std::path::PathBuf;

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
        .get_matches()
}

/// Load a GFA file from the provided path.
pub fn load_gfa<P>(path: P) -> Result<GFA<Vec<u8>, ()>>
where
    P: AsRef<std::path::Path>,
{
    let parser = GFAParser::new();
    let gfa = parser.parse_file(path.as_ref()).with_context(|| {
        format!(
            "Failed to parse GFA from path: {:?}",
            path.as_ref().as_os_str()
        )
    })?;
    Ok(gfa)
}

fn main() -> Result<()> {
    let args = cli();
    let gfa_file = args.get_one::<PathBuf>("GFA").expect("GFA required");
    let gfa = load_gfa(gfa_file).context("Failed to load GFA file")?;

    let repeat_node_size_limit = *args.get_one::<usize>("repeat").unwrap();
    let neighboring_node_minimum = *args.get_one::<usize>("neighbor").unwrap();
    let in_out_threshold = *args.get_one::<usize>("inout").unwrap();

    let segment_sizes: HashMap<Vec<u8>, usize> = gfa
        .segments
        .iter()
        .map(|segment| (segment.name.clone(), segment.sequence.len()))
        .collect();

    let mut edge_map: HashMap<Vec<u8>, HashMap<Orientation, Vec<(Vec<u8>, Orientation)>>> =
        HashMap::new();

    for link in &gfa.links {
        let from = link.from_segment.clone();
        let to = link.to_segment.clone();
        let from_orient = link.from_orient;
        let to_orient = link.to_orient;

        edge_map
            .entry(from.clone())
            .or_default()
            .entry(from_orient)
            .or_default()
            .push((to.clone(), to_orient));

        edge_map
            .entry(to.clone())
            .or_default()
            .entry(to_orient)
            .or_default()
            .push((from.clone(), from_orient));
    }

    let mut repeat_candidates = Vec::new();

    for segment in &gfa.segments {
        let id = segment.name.clone();
        let size = segment.sequence.len();

        // if the segment is too large, skip it
        if size > repeat_node_size_limit {
            continue;
        }

        let neighbor_count: usize = edge_map
            .get(&id)
            .map(|orient_map| orient_map.values().map(|v| v.len()).sum())
            .unwrap_or(0);

        if neighbor_count < in_out_threshold * 2 {
            continue;
        }

        let mut valid_neighbors = true;
        for orient_neighbors in edge_map.get(&id).unwrap_or(&HashMap::new()).values() {
            for (neighbor_id, _) in orient_neighbors {
                if let Some(&neighbor_size) = segment_sizes.get(neighbor_id) {
                    if neighbor_size < neighboring_node_minimum {
                        valid_neighbors = false;
                        break;
                    }
                } else {
                    valid_neighbors = false;
                    break;
                }
            }
            if !valid_neighbors {
                break;
            }
        }

        if valid_neighbors {
            repeat_candidates.push(id);
        }
    }

    if !repeat_candidates.is_empty() {
        println!("ID\tSize");
    }

    for node in repeat_candidates {
        let node_name = std::str::from_utf8(&node)?;
        let segment_size = segment_sizes.get(&node).unwrap();

        println!("{}\t{}", node_name, segment_size);
    }

    Ok(())
}
