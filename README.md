# Scripts to Analyze rDNA from Oxford Nanopore Sequencing Data

Analyze the structure, methylation, variation of rDNA repeats.
IMPORTANT: THROUGHOUT THE PROJECT, QUERY POSITION IS DEFINED IN TERMS OF THE DIRECTION OF THE ORIGINAL BASECALLED READ!

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [Dotplot Visualization](#dotplot-visualization)
  - [R-Repeat Analysis](#r-repeat-analysis)
  - [Methylation Analysis](#methylation-analysis)

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
plot.single_read_plot_structure_fast5(reader, split_analyzed_read, 'dir_name', met=True, pdf=False)
```

#### Using Dorado BAM Files
The most efficient way is to first analyze the data and then plot specific reads:

```python
import dorado_bam_io
import plot

# Analyze BAM files once to get all read data
with dorado_bam_io.BamAnalyzer("mapped.bam", "methylation_called.bam") as analyzer:
    reads_data = analyzer.analyze()

# Plot specific reads using the pre-analyzed data
read_id = "your_read_id_here"
plot.single_read_plot_structure_bam(read_id, reads_data, 'output_dir')
```

#### Batch Plot Multiple Reads
For efficiently plotting multiple reads from the same BAM files:

```python
import plot

# Option 1: Plot specific reads by ID

read_ids = ["read1", "read2", "read3"]
plot.batch_plot_reads_from_bams(
    read_ids, 
    "mapped.bam", 
    "methylation_called.bam", 
    'output_dir',
    met=True  # Include methylation data
)

# Option 2: Plot all reads in BAM files with filtering options
plot.batch_plot_all_reads_from_bams(
    "mapped.bam",
    "methylation_called.bam",
    'output_dir',
    max_reads=100,               # Limit number of reads to process
    min_quality=15,              # Skip reads with average quality below 15
    skip_strand_switch=True,     # Skip reads with strand switching
    pdf=False                    # Save as PNG (default)
)
```

#### Legacy Method (Backward Compatible)
The original function is maintained for backward compatibility:

```python
# For FAST5 files
plot.signle_read_plot_structure(
    f5reader=reader, 
    analyzed_read=split_analyzed_read, 
    savedir='dir_name'
)

# For BAM files
plot.signle_read_plot_structure(
    dorado_bam="methylation_called.bam", 
    mapped_bam="mapped.bam", 
    read_id="your_read_id", 
    savedir='dir_name'
)
```

### R-Repeat Analysis
The r-repeat region contains variable number of repeats. This functionality helps analyze r-repeat composition in reads.

#### Create R-Repeat References
Generate different combinations of r-repeat reference sequences.
```python
from r_repeat_analysis import make_r_repeat_refs
make_r_repeat_refs()
```

#### Map Reads to R-Repeat References
Map reads to r-repeat references and identify the best matching r-repeat pattern for each read.
```python
from r_repeat_analysis import map_r_repeats
map_r_repeats("path/to/reads.fastq", "output_summary_basename")
```

#### Plot R-Repeat Distributions
Visualize the distribution of r-repeat patterns across samples.
```python
from r_repeat_analysis import plot_r_repeat_distributions
plot_r_repeat_distributions("path/to/summary_directory", "output_plot.pdf")
```

### Methylation Analysis
Analyze methylation patterns in reads mapped to rDNA.

#### Using BAM Files with Direct Methylation Analysis
Use the BamAnalyzer to process both a methylation-called BAM and an aligned BAM for the same reads.
```python
import dorado_bam_io

with dorado_bam_io.BamAnalyzer("mapped.bam", "methylation_called.bam") as analyzer:
    reads_data = analyzer.analyze()
    
    # Access read-level data including methylation
    for read_id, read_data in reads_data.items():
        if read_data.has_strand_switch:
            continue  # Skip reads with strand switching
            
        methylation_data = read_data.methylation_data
        # Process methylation data
```

#### Methylation Summary by Reference Position
Create a methylation summary mapped to reference positions, useful for studying site-specific methylation patterns.
```python
from methylation_analysis import make_methylation_summary_in_reference_positions_with_aligned_bam

read_methylation_dict, global_summary = make_methylation_summary_in_reference_positions_with_aligned_bam(
    "aligned.bam", "methylation_called.bam", "reference.fasta"
)

# Access per-read methylation patterns at reference positions
for read_id, data in read_methylation_dict.items():
    for alignment in data['alignments']:
        methylation_at_ref_positions = alignment['methylation']
        # Process methylation data
```

#### Analyze Methylation in Coding Regions
Analyze methylation patterns specifically in coding regions (4-16kb).
```python
from methylation_analysis import methylation_per_coding

# Analyze and create a histogram of methylation scores
scores = methylation_per_coding(
    dorado_bam_path="methylation_called.bam",
    mapped_bam_path="aligned.bam",
    output_file="methylation_histogram.png"
)
```

#### Methylation Binning for Visualization
Process methylation data into bins for easier visualization.
```python
from methylation_analysis import methylation_binning

# Get binned methylation data (x positions and y values)
x, y = methylation_binning(
    window=200,                   # Bin size
    methylation_data=data,        # List of (position, score) tuples
    method='threshold'            # 'average' or 'threshold'
)
```