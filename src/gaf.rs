use anyhow::{Error, Result};
use bstr::{io::*, ByteSlice};
use gfa::{
    gafpaf::{parse_gaf, GAFPath, GAFStep},
    gfa::Orientation,
    optfields::OptField,
};
use std::{
    borrow::BorrowMut,
    collections::{HashMap, HashSet},
    fmt::Display,
    fs::File,
    io::BufReader,
    path::PathBuf,
};

// Example paths from one segment
// define 'repeat' as always on +ve segment
// enumerate In 1, In 2, Out 1, Out 2
// in both forward and reverse orientations
//
// u25	132	<u28<u25>u27 (reverse)
// u25	126	<u27>u25<u28 (forward)
// u25	126	<u27>u25>u28 (forward)
// u25	120	>u28<u25>u26 (reverse)
// u25	119	>u28<u25>u27 (reverse)
// u25	115	<u26>u25<u28 (forward)
// u25	107	<u28<u25>u26 (reverse)
// u25	101	<u26>u25>u28 (forward)

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

    let mut paths: Vec<(&(String, String), &i32)> = paths.iter().collect();
    paths.sort_by(|a, b| b.1.cmp(a.1));

    let paths = Paths::from_vec(paths).split_into_repeats();

    for paths in paths.clone() {
        for ((s, p), c) in paths.paths {
            println!("{p}:{c}");
        }
    }

    output_repeat_lines(paths);

    Ok(())
}

// need a function to:
// - designate 2 in nodes, 2 out nodes
// - identify if a path is reverse of another path
fn output_repeat_lines(all_paths: Vec<Paths>) {
    for paths in all_paths {
        let mut node_checker = Vec::new();
        for ((s, p), c) in paths.paths.clone() {
            for ((s2, p2), c2) in paths.paths.iter().skip(1) {
                let is_reverse = p.is_reverse(p2.clone());
                // if we hit a reverse path, combine them
                // and push to the node checker
                if is_reverse
                    && !node_checker.contains(&p.to_string())
                    && !node_checker.contains(&p2.to_string())
                {
                    println!("{}:{}, {}:{}", p, c, p2, c2);
                    node_checker.push(p.to_string());
                    node_checker.push(p2.to_string());
                }
            }
        }
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

    // function to split the vec into multipe subsets, each of
    // length 8, corresponding to the repeat node
    fn split_into_repeats(&self) -> Vec<Paths> {
        let mut out: Vec<Paths> = Vec::new();
        let mut inner = Paths::new(Vec::new());

        // while the node ID is the same, push to inner
        // then push inner to out and reset inner
        for (i, ((id, path), count)) in self.paths.iter().enumerate() {
            // if we have an empty inner, just push to it
            if inner.paths.is_empty() {
                inner.paths.push(((id.clone(), path.clone()), *count));
            } else {
                let last = inner.paths.last().unwrap();
                if last.0 .0 == id.clone() {
                    inner.paths.push(((id.clone(), path.clone()), *count));
                } else {
                    out.push(inner);
                    inner = Paths::new(Vec::new());
                    inner.paths.push(((id.clone(), path.clone()), *count));
                }
            }
        }

        // if there's only one id present, we never reach to pushing to
        // out, so account for that here
        if !inner.paths.is_empty() {
            out.push(inner);
        }

        out
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
    fn is_reverse(&self, other: Path) -> bool {
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
        assert!(path.is_reverse(path2));
    }

    #[test]
    fn test_is_reverse_2() {
        let s = ">u27>u29<u30";
        let s2 = "<u30<u29<u27";

        let path = string_to_path(s.to_string()).unwrap();
        let path2 = string_to_path(s2.to_string()).unwrap();

        assert!(!path.is_reverse(path2));
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
}
