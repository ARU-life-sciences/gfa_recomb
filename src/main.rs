use anyhow::{Context, Result};
use gfa::gfa::{Link, GFA};
use gfa::optfields::OptFields;
use gfa::parser::GFAParser;
use std::collections::HashMap;

const REPEAT_NODE_SIZE_LIMIT: usize = 10_000;
const NEIGHBORING_NODE_MINIMUM: usize = 5_000;
const IN_OUT_THRESHOLD: usize = 2;

const VERSION: &str = "0.1.0";

/// Given a path, load the GFA into a `GFA` struct.
pub fn load_gfa<T, P>(path: P) -> Result<GFA<Vec<u8>, T>>
where
    T: OptFields,
    P: AsRef<std::path::Path>,
{
    let parser = GFAParser::new();
    let gfa = parser
        .parse_file(path.as_ref())
        .with_context(|| {
            format!(
                "Failed to parse GFA from path: {:?}",
                path.as_ref().as_os_str()
            )
        })
        .context("Failed to parse GFA file")?;
    Ok(gfa)
}

fn main() -> Result<()> {
    // check there is only 1 argument
    if std::env::args().count() != 2 {
        eprintln!("Detect bidirectionally bifurcated repeats in mitochondrial genomes.");
        eprintln!("Paths through the focal node ignore orientation.");
        eprintln!("Version: {}", VERSION);
        eprintln!("Usage: {} <GFA>", std::env::args().next().unwrap());
        std::process::exit(1);
    }

    // Path to the GFA file
    let gfa_file = std::env::args()
        .nth(1)
        .expect("Please provide the path to the GFA file as an argument");

    // Open and parse the GFA file
    let gfa = load_gfa(gfa_file).context("Failed to load GFA file")?;

    // Collect segment sizes into a HashMap
    let segment_sizes: HashMap<Vec<u8>, usize> = gfa
        .segments
        .iter()
        .map(|segment| (segment.name.clone(), segment.sequence.len()))
        .collect();

    // Build adjacency lists
    let mut incoming: HashMap<Vec<u8>, Vec<&Link<Vec<u8>, ()>>> = HashMap::new();
    let mut outgoing: HashMap<Vec<u8>, Vec<&Link<Vec<u8>, ()>>> = HashMap::new();

    for link in &gfa.links {
        let from = link.from_segment.clone();
        let to = link.to_segment.clone();

        // Handle orientation if necessary
        outgoing.entry(from).or_default().push(link);
        incoming.entry(to).or_default().push(link);
    }

    // Now find repeat candidates based on the criteria
    let mut repeat_candidates = Vec::new();

    for segment in &gfa.segments {
        let id = segment.name.clone();
        let size = segment.sequence.len();

        // 1. Check size ≤ 10 kb
        if size > REPEAT_NODE_SIZE_LIMIT {
            // eprintln!(
            //     "Skipping Segment ID: {}, size: {} (too large)",
            //     std::str::from_utf8(&id)?,
            //     size
            // );
            continue;
        }

        // 2. Check 2 incoming and 2 outgoing edges
        // 2. Check 2 unique incoming and 2 unique outgoing edges
        let in_count = incoming
            .get(&id)
            .map(|v| {
                v.iter()
                    .map(|link| link.from_segment.clone())
                    .collect::<Vec<_>>()
                    .len()
            })
            .unwrap_or(0);

        let out_count = outgoing
            .get(&id)
            .map(|v| {
                v.iter()
                    .map(|link| link.to_segment.clone())
                    .collect::<Vec<_>>()
                    .len()
            })
            .unwrap_or(0);

        if in_count + out_count < IN_OUT_THRESHOLD * 2 {
            // eprintln!(
            //     "Skipping Segment ID: {}, size: {} (in: {}, out: {})",
            //     std::str::from_utf8(&id)?,
            //     size,
            //     in_count,
            //     out_count
            // );
            continue;
        }

        // 3. Check all connected nodes are ≥ 10 kb
        let mut valid_neighbors = true;

        let connected_incoming = incoming
            .get(&id)
            .map(|links| {
                links
                    .iter()
                    .map(|l| l.from_segment.clone())
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();

        let connected_outgoing = outgoing
            .get(&id)
            .map(|links| {
                links
                    .iter()
                    .map(|l| l.to_segment.clone())
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();

        for neighbor_id in connected_incoming.iter().chain(connected_outgoing.iter()) {
            // let neighbor_size = segment_sizes.get(neighbor_id);
            // eprintln!(
            //     "  Neighbor ID: {}, size: {}",
            //     std::str::from_utf8(neighbor_id)?,
            //     neighbor_size.unwrap_or(&0)
            // );
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

        if valid_neighbors {
            repeat_candidates.push(id);
        }
    }

    // Output repeat candidates
    for node in repeat_candidates {
        let node_name = std::str::from_utf8(&node)?;
        let segment_size = segment_sizes.get(&node).unwrap();

        let connected_incoming = incoming.get(&node).cloned().unwrap_or(Vec::new());
        let connected_incoming_str = connected_incoming
            .iter()
            .map(|l| std::str::from_utf8(&l.from_segment).unwrap())
            .collect::<Vec<_>>()
            .join(",");

        let connected_outgoing = outgoing.get(&node).cloned().unwrap_or(Vec::new());
        let connected_outgoing_str = connected_outgoing
            .iter()
            .map(|l| std::str::from_utf8(&l.to_segment).unwrap())
            .collect::<Vec<_>>()
            .join(",");

        println!("ID\tSize\tIncoming\tOutgoing");
        println!(
            "{}\t{}\t{}\t{}",
            node_name, segment_size, connected_incoming_str, connected_outgoing_str
        );

        // === ADDED BANDAGE-STYLE PATH ENUMERATOR ===
        println!("Possible paths:");

        let preds: Vec<_> = connected_incoming
            .iter()
            .map(|l| std::str::from_utf8(&l.from_segment).unwrap())
            .collect();

        let succs: Vec<_> = connected_outgoing
            .iter()
            .map(|l| std::str::from_utf8(&l.to_segment).unwrap())
            .collect();

        if !preds.is_empty() && !succs.is_empty() {
            for pred in &preds {
                for succ in &succs {
                    println!("{} -> {} -> {}", pred, node_name, succ);
                }
            }
        } else if !preds.is_empty() {
            for pred in &preds {
                println!("{} -> {}", pred, node_name);
            }
        } else if !succs.is_empty() {
            for succ in &succs {
                println!("{} -> {}", node_name, succ);
            }
        } else {
            println!("{} has no connected paths.", node_name);
        }

        println!();
        // === END BANDAGE-STYLE PATH ENUMERATOR ===
    }
    Ok(())
}
