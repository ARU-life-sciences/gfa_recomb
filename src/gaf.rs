use anyhow::Result;
use bstr::{io::*, ByteSlice};
use gfa::{
    gafpaf::{parse_gaf, GAFPath, GAFStep},
    gfa::Orientation,
    optfields::OptField,
};
use std::f64::consts::LN_2;
use std::io::{self, Write};
use std::{collections::HashMap, fmt::Display, fs::File, io::BufReader, path::PathBuf};

/// Parses a GAF file and extracts all 3-node paths through specified repeat nodes.
/// Counts and groups them by focal repeat segment for recombination analysis.
pub fn count_gaf_paths(gaf_path: PathBuf, nodes: Vec<String>) -> Result<()> {
    let file = File::open(gaf_path).unwrap();
    let lines = BufReader::new(file).byte_lines();

    // store the paths
    let mut paths = HashMap::new();

    for (i, line) in lines.enumerate() {
        let line = line?;
        let fields = line.split_str(b"\t");
        if let Some(gaf) = parse_gaf::<_, Vec<OptField>>(fields) {
            // get the path length
            let path = gaf.path;
            match &path {
                GAFPath::StableId(_) => continue, // don't care about this
                GAFPath::OrientIntv(vec) => {
                    // only interested in those paths of length 3... if
                    // we had longer read data (e.g. Nanopore), we might change this
                    if vec.len() == 3 {
                        let mut node = Default::default();
                        if match vec[1].clone() {
                            GAFStep::SegId(_, id) => {
                                node = String::from_utf8(id.to_vec()).unwrap();
                                nodes.contains(&node)
                            }
                            GAFStep::StableIntv(_, _, _, _) => false,
                        } {
                            *paths.entry((node, path.to_string())).or_insert(0) += 1;
                        } else {
                            continue;
                        }
                    }
                }
            };
        } else {
            eprintln!("Error parsing GAF line {}", i);
            std::process::exit(1);
        }
    }

    let paths: Vec<(&(String, String), &i32)> = paths.iter().collect();

    let paths = Paths::from_vec(paths).split_into_repeats();

    let entropy = compute_path_entropy(&paths);
    output_repeat_lines(paths);

    // format entropy
    let mut stdout = io::stdout();
    let _ = writeln!(stdout, "\nRepeat node\tPath count\tEntropy");
    for (node, path_length, entropy) in entropy.2 {
        let _ = writeln!(stdout, "{}\t{}\t{:.3}", node, path_length, entropy);
    }
    let _ = writeln!(stdout, "\nMean entropy: {:.3}", entropy.0);
    let _ = writeln!(stdout, "Total entropy: {:.3}", entropy.1);

    Ok(())
}

/// Compute the Shannon entropy of path usage for each repeat node.
/// This reflects the diversity of path usage through each focal repeat.
///
/// For each repeat node:
/// - Count all distinct 3-node paths
/// - Sum their coverages
/// - Compute entropy:
///   H = -sum(p_i * log2(p_i)) over all paths i
///
/// Returns:
/// - (mean_entropy, Vec<(repeat_id, path_count, entropy)>)
fn compute_path_entropy(groups: &[Paths]) -> (f64, f64, Vec<(String, usize, f64)>) {
    let mut entropies = Vec::new();

    for group in groups {
        if group.paths.is_empty() {
            continue;
        }

        let repeat_id = &group.paths[0].0 .0; // all paths share this repeat node
        let total_cov: f64 = group.paths.iter().map(|(_, c)| *c as f64).sum();

        if total_cov == 0.0 {
            continue;
        }

        let mut entropy = 0.0;
        for (_, cov) in &group.paths {
            let p = *cov as f64 / total_cov;
            if p > 0.0 {
                entropy -= p * (p.ln() / LN_2);
            }
        }

        entropies.push((repeat_id.clone(), group.paths.len(), entropy));
    }

    let total_entropy = entropies.iter().map(|(_, _, e)| *e).sum::<f64>();
    let mean_entropy = if !entropies.is_empty() {
        entropies.iter().map(|(_, _, e)| *e).sum::<f64>() / entropies.len() as f64
    } else {
        0.0
    };

    (mean_entropy, total_entropy, entropies)
}

/// Compute the Recombination Complexity Index (RCI) from a list of recombination path pairs.
///
/// Each tuple in `revcomps` represents a pair of paths that are reverse complements,
/// along with their coverage values. The RCI summarizes both:
///
/// 1. **Coverage balance** between the forward and reverse paths
/// 2. **Path diversity** through repeat nodes (more paths = higher isomeric potential)
///
/// The RCI is computed as:
///
/// ```text
///     RCI = (1 / R) * Σ [S_r * log2(P_r)]
/// ```
///
/// Where:
/// - `R` is the number of repeat nodes (focal segments)
/// - `P_r` is the number of distinct paths through repeat node `r`
/// - `S_r` is the average recombination score for repeat `r`, where:
///   `recomb_score = 2 * min(rel_cov1, rel_cov2)`
///
/// # Arguments
/// * `revcomps` - A vector of tuples (path1, cov1, path2, cov2)
///
/// # Returns
/// * `f64` - The recombination complexity index (RCI)
///
fn compute_rci(revcomps: &Vec<(String, i32, String, i32)>) -> f64 {
    // Map of repeat node ID to its recombination scores and path count
    let mut repeat_groups: HashMap<String, (Vec<f64>, usize)> = HashMap::new();

    for (p1, cov1, p2, cov2) in revcomps {
        let total_cov = *cov1 as f64 + *cov2 as f64;
        if total_cov == 0.0 {
            continue; // skip zero-coverage cases
        }

        let rel1 = *cov1 as f64 / total_cov;
        let recomb_score = 2.0 * rel1.min(1.0 - rel1);

        // Parse the repeat node from the middle segment of either path
        if let (Some(repeat_id), path_count) = extract_repeat_info(p1, p2) {
            repeat_groups
                .entry(repeat_id)
                .and_modify(|(scores, count)| {
                    scores.push(recomb_score);
                    *count += path_count;
                })
                .or_insert((vec![recomb_score], path_count));
        }
    }

    // Compute average S_r * log2(P_r) for each repeat
    let mut total_score = 0.0;
    let mut n_repeats = 0;

    for (_repeat_id, (scores, path_count)) in repeat_groups {
        if scores.is_empty() || path_count <= 1 {
            continue;
        }

        let s_r: f64 = scores.iter().sum::<f64>() / scores.len() as f64;
        let p_r: f64 = path_count as f64;

        let score = s_r * (p_r.ln() / LN_2); // log2
        total_score += score;
        n_repeats += 1;
    }

    if n_repeats > 0 {
        total_score / n_repeats as f64
    } else {
        0.0
    }
}

/// Extract the focal repeat segment and number of distinct paths from two path strings.
///
/// Assumes paths are of the form `<seg1>seg2<seg3` or similar.
///
/// Returns:
/// - `Some(repeat_id)` and `2` if both paths share the same central repeat segment
/// - `(None, 0)` if paths are malformed or mismatched
fn extract_repeat_info(p1: &str, p2: &str) -> (Option<String>, usize) {
    let repeat_1 = extract_middle_segment(p1);
    let repeat_2 = extract_middle_segment(p2);

    if let (Some(r1), Some(r2)) = (repeat_1, repeat_2) {
        if r1 == r2 {
            return (Some(r1), 2); // this pair contains 2 paths
        }
    }
    (None, 0)
}

/// Extract the middle (repeat) segment from a path string.
/// Assumes the path is made of three segments separated by direction indicators `<` or `>`.
///
/// Example:
/// - Input: `"<u27>u25<u28"` → Output: `Some("u25")`
fn extract_middle_segment(path: &str) -> Option<String> {
    let split: Vec<&str> = path.split(['<', '>']).filter(|s| !s.is_empty()).collect();
    if split.len() == 3 {
        Some(split[1].to_string())
    } else {
        None
    }
}

/// Identifies reverse-complement path pairs within each repeat group, computes their
/// recombination score, prints all path pairs with scores, and reports:
/// - The mean recombination potential
/// - The recombination complexity index (RCI)
fn output_repeat_lines(all_paths: Vec<Paths>) {
    let mut revcomps = Vec::new();
    let mut stdout = io::stdout();

    for paths in all_paths {
        let mut node_checker = Vec::new();

        let path_vec = paths.paths;
        let path_vec: Vec<_> = path_vec.to_vec();

        for ((_, p), c) in path_vec.iter() {
            for ((_, p2), c2) in path_vec.iter().skip(1) {
                let is_reverse = p.is_reverse(p2);
                // if we hit a reverse path, combine them
                // and push to the node checker
                if is_reverse
                    && !node_checker.contains(&p.to_string())
                    && !node_checker.contains(&p2.to_string())
                {
                    revcomps.push((p.to_string(), *c, p2.to_string(), *c2));
                    node_checker.push(p.to_string());
                    node_checker.push(p2.to_string());
                }
            }
        }
    }
    if !revcomps.is_empty() {
        let _ = writeln!(stdout, "path_1\tcov_1\tpath_2\tcov_2\trecomb_score");
        let mut recomb_scores = Vec::new();

        for (p1, cov1, p2, cov2) in &revcomps {
            let total = *cov1 as f64 + *cov2 as f64;
            if total == 0.0 {
                continue; // avoid division by zero
            }
            let rel1 = *cov1 as f64 / total;
            let rel2 = *cov2 as f64 / total;
            let score = 2.0 * rel1.min(rel2);
            recomb_scores.push(score);
            let _ = writeln!(stdout, "{}\t{}\t{}\t{}\t{:.3}", p1, cov1, p2, cov2, score);
        }

        let recomb_potential: f64 = if !recomb_scores.is_empty() {
            recomb_scores.iter().sum::<f64>() / recomb_scores.len() as f64
        } else {
            0.0
        };

        let _ = writeln!(stdout, "\nRecombination potential: {:.3}", recomb_potential);
        let rci = compute_rci(&revcomps);
        let _ = writeln!(stdout, "RCI: {:.3}", rci);
    }
}

#[derive(Debug, Clone)]
struct Paths {
    // node ID, path, count of the path
    paths: Vec<((String, Path), i32)>,
}

impl Paths {
    fn new(inner: Vec<((String, Path), i32)>) -> Self {
        Self { paths: inner }
    }
    fn from_vec(vec: Vec<(&(String, String), &i32)>) -> Self {
        let mut paths = Vec::new();
        for ((id, path), count) in vec {
            let path = string_to_path(path.to_string()).unwrap();
            paths.push(((id.clone(), path), *count));
        }

        Self::new(paths)
    }

    /// Groups path+coverage entries by focal repeat node ID.
    /// Assumes the list is sorted and that paths from the same node appear contiguously.
    fn split_into_repeats(&self) -> Vec<Paths> {
        let mut grouped: HashMap<String, Vec<((String, Path), i32)>> = HashMap::new();

        for entry in &self.paths {
            let repeat_id = entry.0 .0.clone();
            grouped.entry(repeat_id).or_default().push(entry.clone());
        }

        grouped.into_values().map(Paths::new).collect()
    }
}

#[derive(Clone, Debug)]
struct Segment {
    orientation: Orientation,
    segid: String,
}

#[derive(Debug, Clone)]
struct Path([Segment; 3]);

impl Display for Path {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut out = String::new();
        for segment in self.0.iter() {
            out.push_str(&format!(
                "{}{}",
                match segment.orientation {
                    Orientation::Forward => ">",
                    Orientation::Backward => "<",
                },
                segment.segid
            ));
        }
        write!(f, "{}", out)
    }
}

impl Path {
    fn from(&self) -> Segment {
        self.0[0].clone()
    }

    fn repeat(&self) -> Segment {
        self.0[1].clone()
    }

    fn to(&self) -> Segment {
        self.0[2].clone()
    }
    // u25	132	<u28<u25>u27 (reverse)
    // u25	126	<u27>u25>u28 (forward)
    // these two paths are reverses of one another
    fn is_reverse(&self, other: &Path) -> bool {
        // the first and third orientations should differ, but should both have
        // the same segid
        let first_third_o = self.from().orientation != other.to().orientation;
        let first_third_p = self.from().segid == other.to().segid;

        // but also need to compare the third and the first
        let third_first_o = self.to().orientation != other.from().orientation;
        let third_first_p = self.to().segid == other.from().segid;

        // the repeat should be the same, but the orientation should be different
        let repeat_same = self.repeat().segid == other.repeat().segid;
        let repeat_diff = self.repeat().orientation != other.repeat().orientation;

        // return the bool
        first_third_o
            && first_third_p
            && third_first_o
            && third_first_p
            && repeat_same
            && repeat_diff
    }
}

impl Default for Segment {
    fn default() -> Self {
        Segment {
            orientation: Orientation::Forward,
            segid: String::new(),
        }
    }
}

/// Parses a string representation of a 3-segment path (e.g. ">u28<u25>u26") into a Path object.
/// Validates directionality (< or >) and ensures exactly 3 segments. Returns an error if malformed.
fn string_to_path(s: String) -> Result<Path> {
    let mut out: [Segment; 3] = core::array::from_fn(|_| Segment::default());

    // split the string on < or >
    let split: Vec<&str> = s.split(['<', '>']).collect();
    // remove empty strings in vec
    let split: Vec<&str> = split.into_iter().filter(|&x| !x.is_empty()).collect();
    // also collect a vec of the delimiters
    let delimiters: Vec<&str> = s.split(|c| c != '<' && c != '>').collect();
    // remove empty strings in vec
    let delimiters: Vec<&str> = delimiters.into_iter().filter(|&x| !x.is_empty()).collect();

    let mut inner_segment = Segment::default();
    for (i, (path, orientation)) in split.iter().zip(delimiters).enumerate() {
        inner_segment.orientation = match orientation {
            "<" => Orientation::Backward,
            ">" => Orientation::Forward,
            _ => return Err(anyhow::anyhow!("no orientation")), // default
        };
        inner_segment.segid = path.to_string();
        out[i] = inner_segment;

        // reset the inner segment
        inner_segment = Segment::default();
    }

    Ok(Path(out))
}

// tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_parse() {
        let s = "<u28<u25>u27";
        let path = string_to_path(s.to_string()).unwrap();
        assert_eq!(path.from().segid, "u28");
        assert_eq!(path.from().orientation, Orientation::Backward);
        assert_eq!(path.repeat().segid, "u25");
        assert_eq!(path.repeat().orientation, Orientation::Backward);
        assert_eq!(path.to().segid, "u27");
        assert_eq!(path.to().orientation, Orientation::Forward);
    }

    #[test]
    fn test_is_reverse() {
        let s = "<u28<u25>u27";
        let path = string_to_path(s.to_string()).unwrap();
        let s2 = "<u27>u25>u28";
        let path2 = string_to_path(s2.to_string()).unwrap();
        assert!(path.is_reverse(&path2));
    }

    #[test]
    fn test_is_reverse_2() {
        let s = ">u27>u29<u30";
        let s2 = "<u30<u29<u27";

        let path = string_to_path(s.to_string()).unwrap();
        let path2 = string_to_path(s2.to_string()).unwrap();

        assert!(!path.is_reverse(&path2));
    }

    // test split into repeats
    // u25	132	<u28<u25>u27 (reverse)
    // u25	126	<u27>u25<u28 (forward)
    // u25	126	<u27>u25>u28 (forward)
    // u25	120	>u28<u25>u26 (reverse)
    // u25	119	>u28<u25>u27 (reverse)
    // u25	115	<u26>u25<u28 (forward)
    // u25	107	<u28<u25>u26 (reverse)
    // u25	101	<u26>u25>u28 (forward)
    #[test]
    fn test_split_repeats() {
        let paths = Paths::from_vec(vec![
            (&("u25".into(), "<u28<u25>u27".into()), &132),
            (&("u25".into(), "<u27>u25<u28".into()), &126),
            (&("u25".into(), "<u27>u25>u28".into()), &126),
            (&("u25".into(), ">u28<u25>u26".into()), &120),
            (&("u25".into(), ">u28<u25>u27".into()), &119),
            (&("u25".into(), "<u26>u25<u28".into()), &115),
            (&("u25".into(), "<u28<u25>u26".into()), &107),
            (&("u25".into(), "<u26>u25>u28".into()), &101),
        ]);

        let split = paths.split_into_repeats();

        assert!(split.len() == 1);
    }

    #[test]
    fn test_split_repeats_multiple() {
        // +u30-u29-u26 -
        // -u30-u29-u26 -
        // +u26+u29-u30 -
        // -u30-u29-u27 -
        // +u27+u29+u30 -
        // +u30-u29-u27 -
        // +u26+u29+u30 -
        let paths = Paths::from_vec(vec![
            (&("u25".into(), "<u28<u25>u27".into()), &132),
            (&("u25".into(), "<u27>u25<u28".into()), &126),
            (&("u25".into(), "<u27>u25>u28".into()), &126),
            (&("u25".into(), ">u28<u25>u26".into()), &120),
            (&("u25".into(), ">u28<u25>u27".into()), &119),
            (&("u25".into(), "<u26>u25<u28".into()), &115),
            (&("u25".into(), "<u28<u25>u26".into()), &107),
            (&("u25".into(), "<u26>u25>u28".into()), &101),
            (&("u29".into(), ">u27>u29<u30".into()), &87),
            (&("u29".into(), ">u30<u29<u26".into()), &67),
            (&("u29".into(), "<u30<u29<u26".into()), &66),
            (&("u29".into(), ">u26>u29<u30".into()), &65),
            (&("u29".into(), "<u30<u29<u27".into()), &63),
            (&("u29".into(), ">u27>u29>u30".into()), &61),
            (&("u29".into(), ">u30<u29<u27".into()), &57),
            (&("u29".into(), ">u26>u29>u30".into()), &55),
        ]);

        let split = paths.split_into_repeats();

        assert!(split.len() == 2);
        assert_eq!(split[0].paths.len(), 8);
        assert_eq!(split[1].paths.len(), 8);
    }

    #[test]
    fn test_extract_middle_segment_valid() {
        let p = "<u27>u25<u28";
        assert_eq!(extract_middle_segment(p), Some("u25".to_string()));
    }

    #[test]
    fn test_extract_middle_segment_invalid() {
        let p = ">u27<u25"; // Only two segments
        assert_eq!(extract_middle_segment(p), None);
    }

    #[test]
    fn test_extract_repeat_info_valid_pair() {
        let (id, count) = extract_repeat_info("<u27>u25<u28", ">u26<u25>u27");
        assert_eq!(id, Some("u25".to_string()));
        assert_eq!(count, 2);
    }

    #[test]
    fn test_extract_repeat_info_mismatch() {
        let (id, count) = extract_repeat_info("<u27>u25<u28", ">u26<u24>u27");
        assert_eq!(id, None);
        assert_eq!(count, 0);
    }

    #[test]
    fn test_compute_rci_balanced() {
        let revcomps = vec![
            // repeat node: u66
            (
                "<u67<u66>u65".to_string(),
                100,
                "<u65>u66>u67".to_string(),
                100,
            ),
            // repeat node: u69
            (
                ">u65<u69<u68".to_string(),
                80,
                ">u68>u69<u65".to_string(),
                80,
            ),
        ];

        let rci = compute_rci(&revcomps);
        eprintln!("RCI: {rci}");
        assert!(rci > 0.0);
        assert!((rci - 1.0).abs() < 0.1);

        let revcomps = vec![
            // Repeat node u66 (2 pairs → 4 paths total)
            (
                "<u67<u66>u65".to_string(),
                100,
                "<u65>u66>u67".to_string(),
                100,
            ),
            (
                "<u68<u66>u64".to_string(),
                90,
                "<u64>u66>u68".to_string(),
                90,
            ),
            // Repeat node u69 (2 pairs → 4 paths total)
            (
                ">u65<u69<u68".to_string(),
                80,
                ">u68>u69<u65".to_string(),
                80,
            ),
            (
                ">u64<u69<u67".to_string(),
                70,
                ">u67>u69<u64".to_string(),
                70,
            ),
        ];

        let rci = compute_rci(&revcomps);
        eprintln!("RCI: {rci}");
        assert!(rci > 0.0);
        assert!((rci - 2.0).abs() < 0.1); // 2 repeat nodes × log2(2) × score ~ 1.0
    }

    #[test]
    fn test_compute_rci_unbalanced() {
        let revcomps = vec![
            (
                "<u27>u25<u28".to_string(),
                190,
                ">u26<u25>u27".to_string(),
                10,
            ), // very unbalanced
        ];

        let rci = compute_rci(&revcomps);
        assert!(rci < 0.5);
    }

    #[test]
    fn test_compute_rci_empty() {
        let revcomps = vec![];
        let rci = compute_rci(&revcomps);
        assert_eq!(rci, 0.0);
    }

    #[test]
    fn test_compute_path_entropy_single_node() {
        let paths = Paths::new(vec![
            (
                ("u66".into(), string_to_path("<u67<u66>u65".into()).unwrap()),
                100,
            ),
            (
                ("u66".into(), string_to_path("<u65>u66>u67".into()).unwrap()),
                100,
            ),
            (
                ("u66".into(), string_to_path("<u68<u66>u64".into()).unwrap()),
                50,
            ),
        ]);

        let (mean_entropy, _total_entropy, details) = compute_path_entropy(&[paths]);

        assert!((mean_entropy > 1.0) && (mean_entropy < 2.0));
        assert_eq!(details.len(), 1);
        assert_eq!(details[0].0, "u66");
        assert_eq!(details[0].1, 3);
        assert!(details[0].2 > 1.0); // entropy should be positive
    }

    #[test]
    fn test_split_into_repeats_merges_by_repeat_id() {
        let paths = Paths::new(vec![
            (
                ("u66".into(), string_to_path("<u67<u66>u65".into()).unwrap()),
                100,
            ),
            (
                ("u69".into(), string_to_path(">u65<u69<u67".into()).unwrap()),
                175,
            ),
            (
                ("u66".into(), string_to_path("<u65>u66>u67".into()).unwrap()),
                180,
            ),
            (
                ("u69".into(), string_to_path(">u67>u69<u65".into()).unwrap()),
                171,
            ),
        ]);

        let merged = paths.split_into_repeats();

        assert_eq!(merged.len(), 2);
        let mut ids: Vec<String> = merged.iter().map(|p| p.paths[0].0 .0.clone()).collect();
        ids.sort();
        assert_eq!(ids, vec!["u66", "u69"]);
    }
}
