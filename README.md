# Scripts to analyze rDNA from Oxford Nanopore seq data

Analyze structure, methylation, variation of rDNA repeat.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)

## Installation

Just clone!

## Usage

### Dotplot visualization
fast5
```python
reader = fast5_loader.SingleFast5Reader("some.fast5", methylation='cpg')
analyzer = seq_analysis.ReadAnalyzer(reader, config.RDNA_REF_HUMAN)
split_analyzed_read = analyzer.split_alignment()
plot.signle_read_plot_structure(reader, split_analyzed_read, 'temp_figs_dir/', met=True, pdf=False)
```
Dorado bam
```python
with BamToSingleReadReader(bam_path) as bam_reader:
    for reader in bam_reader:
        analyzer = seq_analysis.ReadAnalyzer(reader, config.RDNA_REF_HUMAN)
        split_analyzed_read = analyzer.split_alignment()
        plot.signle_read_plot_structure(reader, split_analyzed_read, 'temp_figs_dir/', met=True, pdf=False)
```

### Dorado methylation analysis
First, prepare a .bam called by Dorado and a .bam mapped to rDNA with minimap2 (after converting Dorado output .bam to .fastq file by samtools fastq).