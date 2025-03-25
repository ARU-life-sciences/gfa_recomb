use anyhow::Result;
use bstr::{io::*, ByteSlice};
use gfa::{
    gafpaf::{parse_gaf, GAFPath, GAFStep},
    optfields::OptField,
};
use std::{collections::HashMap, fs::File, io::BufReader, path::PathBuf};

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
                        // remove the first and last elements of the vec
                        let mut check_inner = vec.clone();
                        let _ = check_inner.pop();
                        let _ = check_inner.remove(0);
                        // check if the inners contain any of the nodes
                        let mut node = Default::default();
                        if check_inner.iter().any(|x| match x {
                            GAFStep::SegId(_, id) => {
                                node = String::from_utf8(id.to_vec()).unwrap();
                                nodes.contains(&node)
                            }
                            GAFStep::StableIntv(_, _, _, _) => false,
                        }) {
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

    // sorted values
    let mut paths: Vec<(&(String, String), &i32)> = paths.iter().collect();
    paths.sort_by(|a, b| b.1.cmp(a.1));

    // print headers
    println!("ID\tCount\tPath");
    // print the paths
    for ((id, path), count) in paths {
        println!("{}\t{}\t{}", id, count, path);
    }

    Ok(())
}
