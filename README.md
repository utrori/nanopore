# Scripts to Analyze rDNA from Oxford Nanopore Sequencing Data

Analyze the structure, methylation, variation of rDNA repeats.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)

## Installation

Just clone! Needs numpy, matplotlib and pysam. Needs minimap2 and bwa in system PATH.

## Usage

### Dotplot Visualization
File name is read id by default.
#### Single Read FAST5 Files
Maybe you need to specify basecall number by bc_n.
```python
reader = fast5_loader.SingleFast5Reader("some.fast5", methylation='cpg')
analyzer = seq_analysis.ReadAnalyzer(reader, config.RDNA_REF_HUMAN)
split_analyzed_read = analyzer.split_alignment()
plot.signle_read_plot_structure(reader, split_analyzed_read, 'dir_name', met=True, pdf=False)
```
#### Dorado-Called Multi-Read BAM Files
The BAM file should be Dorado called with methylation setting.
```python
with BamToSingleReadReader(bam_path) as bam_reader:
    for reader in bam_reader:
        analyzer = seq_analysis.ReadAnalyzer(reader, config.RDNA_REF_HUMAN)
        split_analyzed_read = analyzer.split_alignment()
        plot.signle_read_plot_structure(reader, split_analyzed_read, 'dir_name', met=True, pdf=False)
```

### Dorado methylation analysis
Prepare two BAM files: one with Dorado basecalls and methylation calls, and another with the same reads mapped to your rDNA reference using Minimap2 (after converting the Dorado BAM to FASTQ with `samtools fastq`).