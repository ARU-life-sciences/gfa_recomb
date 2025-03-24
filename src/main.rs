//! Bidirected Path Enumerator with Flexible Orientation Traversal
//! This program identifies potential repeat nodes in a GFA graph.

use anyhow::{Context, Result};
use gfa::gfa::{Orientation, GFA};
use gfa::parser::GFAParser;
use std::collections::HashMap;

const REPEAT_NODE_SIZE_LIMIT: usize = 10_000;
const NEIGHBORING_NODE_MINIMUM: usize = 5_000;
const IN_OUT_THRESHOLD: usize = 2;

const VERSION: &str = "0.5.0";

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
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 || args.len() > 3 {
        eprintln!("Bidirected Repeat Path Enumerator");
        eprintln!("Version: {}", VERSION);
        eprintln!("Usage: {} <GFA>", args[0]);
        eprintln!(
            "\nBy default, traversal allows entry and exit through the same orientation.\nUse --strict to enforce bidirected traversal (enter one end, exit the other).\n"
        );
        std::process::exit(1);
    }

    let gfa_file = &args[1];

    if gfa_file == "-h" || gfa_file == "--help" {
        eprintln!("This program identifies potential repeat nodes in a GFA graph.");
        std::process::exit(1);
    }

    let gfa = load_gfa(gfa_file).context("Failed to load GFA file")?;

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

        if size > REPEAT_NODE_SIZE_LIMIT {
            continue;
        }

        let neighbor_count: usize = edge_map
            .get(&id)
            .map(|orient_map| orient_map.values().map(|v| v.len()).sum())
            .unwrap_or(0);

        if neighbor_count < IN_OUT_THRESHOLD * 2 {
            continue;
        }

        let mut valid_neighbors = true;
        for orient_neighbors in edge_map.get(&id).unwrap_or(&HashMap::new()).values() {
            for (neighbor_id, _) in orient_neighbors {
                if let Some(&neighbor_size) = segment_sizes.get(neighbor_id) {
                    if neighbor_size < NEIGHBORING_NODE_MINIMUM {
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

    println!("ID\tSize");

    for node in repeat_candidates {
        let node_name = std::str::from_utf8(&node)?;
        let segment_size = segment_sizes.get(&node).unwrap();

        println!("{}\t{}", node_name, segment_size);
    }

    Ok(())
}
