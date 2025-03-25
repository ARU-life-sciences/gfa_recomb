# Detect potential recombination points in plant mitochondrial genomes

Using only information from incoming and outgoing nodes (ignoring orientation), and with some constraints, we can detect all bidirectionally bifurcating segments in a graph.

```bash
# Simple CLI
gfa_recomb <GFA>
```
## GraphAligner output

Including the `--gaf <GAF>` option iterates over the GAF to find alignments which span a focal node (only paths of length 3 considered at the moment). Example output is below.

```
ID	Count	Path
u75	198	<u76>u75<u71
u75	197	>u71<u75>u76
u75	179	>u73>u75<u71
u75	178	<u76>u75<u76
u75	163	>u76<u75>u76
u75	162	>u76<u75<u73
u75	159	>u71<u75<u73
u75	133	>u73>u75<u76
u77	114	>u72>u77>u78
u77	100	<u78<u77<u72
u77	89	>u78<u77<u72
u77	76	>u72>u77<u78
u77	74	>u74>u77>u78
u77	67	<u78<u77<u74
u77	60	>u78<u77<u74
u77	54	>u74>u77<u78
u70	19	<u71>u70>u73
u70	18	<u73<u70>u71
u70	14	<u71>u70>u74
u70	13	<u74<u70>u71
u70	12	<u74<u70>u72
u70	11	<u73<u70>u72
u70	6	<u72>u70>u73
u70	5	<u72>u70>u74
```

## Recombination metric of the GFA

Using the above info, we can begin to assess the recombination potential of each assembly.
