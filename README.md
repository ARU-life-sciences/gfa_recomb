# Detect potential recombination points in plant mitochondrial genomes

Using only information from incoming and outgoing nodes (ignoring orientation), and with some constraints, we can detect all bidirectionally bifurcating segments in a graph.

```bash
# Simple CLI
gfa_recomb <GFA>
```
## GraphAligner output

Including the `--gaf <GAF>` option iterates over the GAF to find alignments which span a focal node (only paths of length 3 considered at the moment). Example output is below.

An example:

```
gfa_recomb --gaf ./data/Arabidopsis_thaliana.gaf ./data/Arabidopsis_thaliana.mito.gfa
```
Should give the following output.

```
path_1  cov_1   path_2  cov_2   recomb_score
<u67<u66>u65    192     <u65>u66>u67    180     0.968
<u68<u66>u64    168     <u64>u66>u68    162     0.982
>u64<u69<u67    160     >u67>u69<u64    126     0.881
>u65<u69<u68    159     >u68>u69<u65    147     0.961
>u68>u69<u64    152     >u64<u69<u68    123     0.895

Recombination potential: 0.937
RCI: 2.154
```

`path_1` and `path_2` are opposite traversals through the same putative repeat node, with approximately similar coverages in the case of `<u67<u66>u65`.

## Recombination metric of the GFA

I propose a new metric, RCI (recombination complexity index).

This quantifies the potential for repeat-mediated structural variation in a genome assembly. It captures both the abundance and diversity of alternative paths through repeat nodes in a GFA graph.

For a set of nodes, R:

RCI = (1 / |R|) * sum over r [ S_r * log2(P_r) ]

Where:
- |R| is the number of distinct repeat nodes (focal segments)
- P_r is the number of distinct paths through repeat node `r`
- S_r is the average recombination score at repeat node `r`

The recombination score at a repeat node, `r`.
S = 2 * min(cov1 / (cov1 + cov2), cov2 / (cov1 + cov2))

This value ranges from:
- 1.0 → perfectly balanced recombination (equal path support)
- 0.0 → only one path is supported (no recombination signal)
