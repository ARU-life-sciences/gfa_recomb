# Detect potential recombination points in plant mitochondrial genomes

Using only information from incoming and outgoing nodes (ignoring orientation), and with some constraints, we can detect all bidirectionally bifurcating segments in a graph.

```bash
# Simple CLI
gfa_recomb <GFA>
```
## GraphAligner output

Optionally take GraphAligner output (GAF format), which will:

- 
- look for paths of length 3 (e.g. >u2<u3<u7)
 
