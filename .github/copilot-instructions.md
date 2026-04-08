# AGFusion

AGFusion is a Python bioinformatics package for annotating gene fusions from human or mouse genomes. Given two fusion gene partners and their junction coordinates, it predicts the functional effect (in-frame, out-of-frame, etc.), generates FASTA sequences, and produces protein domain/exon architecture visualizations.

## Stack

- **Language**: Python 3.9 - 3.13
- **Formatter**: black (line-length = 100)
- **Key dependencies**: pyensembl, biopython, matplotlib, pandas, numpy<2
- **CLI**: Click-based (`agfusion annotate`, `agfusion batch`, `agfusion download`)
- **Database**: SQLite via `AGFusionDB`, storing Pfam/TMHMM protein domain annotations
- **Tests**: unittest, located in `test/`

## Project Layout

```
agfusion/
  model.py      # Core domain classes: _Gene, Fusion
  database.py   # AGFusionDB — SQLite protein domain store
  parsers.py    # Parsers for fusion-finding algorithm outputs
  plot.py       # Matplotlib visualizations
  cli.py        # Click CLI entry points
  utils.py      # Shared utilities
  exceptions.py # Custom exceptions
test/           # unittest test cases + fixture data
```

## Domain Terminology

- Gene partners: `gene5prime`, `gene3prime`
- Junction coordinates: `junction5prime`, `junction3prime` (genomic positions)
- Supported builds: `hg38`, `hg19`, `mm10`
- Ensembl annotation via `pyensembl.EnsemblRelease`
- Default behavior: canonical isoforms only; `noncanonical=True` enables all
