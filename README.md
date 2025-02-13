# Scripts to Analyze rDNA from Oxford Nanopore Sequencing Data

Analyze the structure, methylation, variation of rDNA repeats.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)

## Installation

Just clone! Needs numpy, matplotlib and pysam. Needs minimap2 and bwa in PATH.

## Usage

### Dotplot Visualization
#### Single Read FAST5 FILE
```python
reader = fast5_loader.SingleFast5Reader("some.fast5", methylation='cpg')
analyzer = seq_analysis.ReadAnalyzer(reader, config.RDNA_REF_HUMAN)
split_analyzed_read = analyzer.split_alignment()
plot.signle_read_plot_structure(reader, split_analyzed_read, 'temp_figs_dir/', met=True, pdf=False)
```
#### Dorado-Called Multi-Read BAM Files
```python
with BamToSingleReadReader(bam_path) as bam_reader:
    for reader in bam_reader:
        analyzer = seq_analysis.ReadAnalyzer(reader, config.RDNA_REF_HUMAN)
        split_analyzed_read = analyzer.split_alignment()
        plot.signle_read_plot_structure(reader, split_analyzed_read, 'temp_figs_dir/', met=True, pdf=False)
```

### Dorado methylation analysis
First, prepare a .bam called by Dorado and a .bam mapped to rDNA with minimap2 (after converting Dorado output .bam to .fastq file by samtools fastq).