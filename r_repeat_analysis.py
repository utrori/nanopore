from pathlib import Path
import seaborn as sns
import pandas as pd
import numpy as np
import os
import subprocess
from natsort import natsorted
import collections
import pysam
import utilities
import io
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def make_r_repeat_refs():
    """
    Create FASTA files for different combinations of r-repeat sequences from human rDNA.
    
    This function extracts r-repeat segments from a human rDNA reference sequence and creates
    multiple FASTA files containing different combinations and copy numbers of these repeats.
    
    The function generates:
    - Files with 4 copies of r-repeats (various combinations)
    - Files with 3 copies of r-repeats (various combinations)
    - Files with 2 copies of r-repeats
    - Files with 1 copy of each r-repeat
    
    Each generated sequence includes 500bp flanking regions (seq_5 and seq_3) on either side
    of the repeat combinations.
    
    The r-repeat regions are defined as:
    - repeat1: positions 13493-14288 in the reference
    - repeat2: positions 14289-15022 in the reference
    - repeat3: positions 15023-15584 in the reference
    
    All output files are saved in the "references/r_repeats/" directory.
    
    Returns:
        None
    """
    human_rDNA_ref= Path("references/human_rDNA.fa")
    seq = ""
    with open(human_rDNA_ref, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq += line.strip()
    repeat1_coordinates = (13493, 14288)
    repeat1_seq = seq[repeat1_coordinates[0]:repeat1_coordinates[1]]
    repeat2_coordinates = (14289, 15022)
    repeat2_seq = seq[repeat2_coordinates[0]:repeat2_coordinates[1]]
    repeat3_coordinates = (15023, 15584)
    repeat3_seq = seq[repeat3_coordinates[0]:repeat3_coordinates[1]]
    seq_5 = seq[12993:13493]
    seq_3 = seq[15584:16084]
    with open("references/r_repeats/r_repeat_4copies_1123.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat1_seq}{repeat1_seq}{repeat2_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_4copies_1223.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat1_seq}{repeat2_seq}{repeat2_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_4copies_1233.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat1_seq}{repeat2_seq}{repeat3_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_3copies_123.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat1_seq}{repeat2_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_3copies_122.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat1_seq}{repeat2_seq}{repeat2_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_3copies_133.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat1_seq}{repeat3_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_3copies_223.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat2_seq}{repeat2_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_3copies_233.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat2_seq}{repeat3_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_3copies_113.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat1_seq}{repeat1_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_2copies_12.fa", 'w') as f:
        f.write(f">r_repeat_2copies\n{seq_5}{repeat1_seq}{repeat2_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_2copies_13.fa", 'w') as f:
        f.write(f">r_repeat_2copies\n{seq_5}{repeat1_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_2copies_23.fa", 'w') as f:
        f.write(f">r_repeat_2copies\n{seq_5}{repeat2_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_1copy_1.fa", 'w') as f:
        f.write(f">r_repeat_2copy\n{seq_5}{repeat1_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_1copy_2.fa", 'w') as f:
        f.write(f">r_repeat_2copy\n{seq_5}{repeat2_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_1copy_3.fa", 'w') as f:
        f.write(f">r_repeat_2copy\n{seq_5}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_0copy.fa", 'w') as f:
        f.write(f">r_repeat_2copy\n{seq_5}{seq_3}\n")

def map_r_repeats(fq, summary_basename):
    """
    Maps reads from a FASTQ file to r-repeat reference sequences and identifies the best matching r-repeat for each read.
    This function:
    1. Maps reads from the input FASTQ file to each r-repeat reference sequence using minimap2
    2. Filters alignments based on read length (>30kb) and alignment score (>500)
    3. For each read, identifies the r-repeat with the highest alignment score
    4. Outputs two summary files:
       - A file listing each read ID and its best matching r-repeat type
       - A distribution file showing counts of reads matching each r-repeat type
    Parameters:
    -----------
    fq : str
        Path to the input FASTQ file containing nanopore reads
    summary_basename : str
        Base filename for the output summary files. Will generate:
        - {summary_basename}_read_id2r_repeat_type.txt: tab-delimited file with read_id and best r-repeat match
        - {summary_basename}_repeat_summary.txt: summary counts of reads per r-repeat type
    Returns:
    --------
    None
        Results are written to output files
    Notes:
    ------
    - Expects r-repeat reference sequences in 'references/r_repeats/' directory
    - Uses utility functions for minimap2 mapping
    - Requires reads to be >30kb with alignment score >500 to be considered
    """

    refs = Path('references/r_repeats').glob('*')
    summary = {} # key: read_id, value: dict with keys as r_repeat ref names and values as the alignment score
    r_repeats = []
    for ref in refs:
        r_repeat = ref.stem
        r_repeats.append(r_repeat)
        
        # Use the original approach with the standard minimap2_mapping function
        sam_path = utilities.minimap2_mapping(fq, ref)
        try:
            with pysam.AlignmentFile(sam_path, "r") as f:
                for read in f:
                    if not read.is_mapped or read.query_length < 30000:
                        continue
                    alignment_score = read.get_tag("AS")
                    if alignment_score < 500:
                        continue
                    if read.query_name not in summary:
                        summary[read.query_name] = {}
                    summary[read.query_name][r_repeat] = alignment_score
        finally:
            # Ensure we clean up the temporary file
            if Path(sam_path).exists():
                os.unlink(sam_path)
    
    each_read_r_repeat_ret = ''
    read_id2max_r_repeat = {}
    for read_id, scores in natsorted(summary.items()):
        max_r_repeat = max(scores, key=lambda k: scores[k])
        #print(read_id, max_r_repeat)
        read_id2max_r_repeat[read_id] = max_r_repeat
        each_read_r_repeat_ret += f'{read_id}\t{max_r_repeat}\n'
    with open(f'{summary_basename}_read_id2r_repeat_type.txt', 'w') as fw:
        fw.write(each_read_r_repeat_ret)
    max_dist = collections.defaultdict(int) 
    for r_repeat in read_id2max_r_repeat.values():
        max_dist[r_repeat] += 1
    summary_ret = ''
    for r_repeat in sorted(max_dist):
        summary_ret += f'{r_repeat}\t{max_dist[r_repeat]}\n'
    with open(f'{summary_basename}_repeat_summary.txt', 'w') as fw:
        fw.write(summary_ret)


def get_each_r_repeat_from_seq():
    fq = "test_files/dorado_output_PSCA0047/PSCA0047_dorado.fastq"
    #fq = "test_files/dorado_output_PSCA0047/"
    r_ex5 = "GCGGCCACCCGGGGTCCCGGCCCTCGCGCGTCCTTCCTCCTCGCTCCTCCGCACGGGTCGACCAGCAGACCGCGGGTGGTGGGCGGCGGGCGGCGAGGCCCCACGGGGCGTCCGCGCACCCGGCCGACCTCCGCTCGTGACCTCTCCTCGGTCGGGCCTCCGGGGTCGACCGCCTGCCGCCCGCGGGCG"
    r_ex3 = "GCTGCTGCTGCTGCCTCTGCCTCCACGGTTCAAGCAAACAGCAAGTTTTCTATTTCGAGTAAAGACGTAATTTCACCATTTTGGCCGGGCTGGTCTCGAACTCCCGACCTAGTGATCCGCCCGCCTCGGCCTCCCAAAGACTGCTGGGAGTACAGATGTGAGCCACCATGCCCGGCCGATTCCTTCCTTTTTTCAATCTTATTTTCTGAACGCTGCCGTGTATGAACATACATCTACACATAC"
    consensus_5 = 'TGAGACTCAGCCGGCGTCTCGCCGTGTCCCGGGTCGACCGGCGGGCCTTCTCCACCGAGCGGCGTGTAGGAGTGCCCGTCGGGACGAACCGCAACCGGAGCGTCCCCGTCTCGGTCGGCACCTCCGGGGTCGACCAGCTGCCGCCCGCGAGCTCCGGACTTAGCCGGCGCCTGCACGTGTCCCGGGTCGACCAGCAGGCGGCCGCCGGACGCTGCGGCGCACCGACGCGAGGGCGTCGATTCCCGTTCGCGCGCCCGCGACCTCCACCGGCCTCGGCCCGCGGTGGAGCTGGGACCACGCGGAACTCCCTCTCTCACATTTTTTTCAGCCCCACCGCGAGTTTGCGTCCGCGGGACTTTTAAGAGGGAGTCACTGCTGCCGTCAGCCAGTAATGCTTCCTCCTTTTTTGCTTTT'
    consensus_3 = 'TCCTTGGTGCCTTCTCGGCTC'
    r_ex5_file = 'references/r_repeat_ex5.fa'
    r_ex3_file = 'references/r_repeat_ex3.fa'
    consensus_5_file = 'references/r_repeat_5_consensus.fa'
    consensus_3_file = 'references/r_repeat_3_consensus.fa'
    with open(consensus_5_file, 'w') as fw:
        fw.write(f'>r_repeat_5_consensus\n{consensus_5}')
    with open(consensus_3_file, 'w') as fw:
        fw.write(f'>r_repeat_3_consensus\n{consensus_3}')
    with open(r_ex5_file, 'w') as fw:
        fw.write(f'>r_repeat_ex5\n{r_ex5}')
    with open(r_ex3_file, 'w') as fw:
        fw.write(f'>r_repeat_ex3\n{r_ex3}')
    
    # Use the original minimap2_mapping function
    sam_5 = utilities.minimap2_mapping(fq, consensus_5_file)
    sam2 = utilities.bwa_mapping(fq, consensus_3_file)  
    sam_ex5 = utilities.minimap2_mapping(fq, r_ex5_file)
    sam_ex3 = utilities.minimap2_mapping(fq, r_ex3_file)
    
    read_id2alignments = collections.defaultdict(list)
    
    # Process the alignment using the standard approach
    with pysam.AlignmentFile(sam_ex3, "r") as f:
        for read in f:
            if read.is_mapped:
                read_id2alignments[read.query_name].append((read.reference_start, read.reference_end))
    
    # Clean up temporary files
    for sam_file in [sam_5, sam2, sam_ex5, sam_ex3]:
        if Path(sam_file).exists():
            os.unlink(sam_file)
    
    for rid, alignments in read_id2alignments.items():
        print(alignments)


def analyze_hpgp_r_repeats():
    dorado_bams = Path("/media/owner/bbeedd59-668c-4d67-a65b-442d3272c3bd/hpgp/dorado_bc/").glob("*.bam")
    for dorado_bam in dorado_bams:
        sample_name = dorado_bam.stem
        print(sample_name)
        temp_fastq = Path("temp_files") / f"{sample_name}_temp.fastq"
        subprocess.run(f"samtools fastq -@ 10 {dorado_bam} > {temp_fastq}", shell=True)
        map_r_repeats(temp_fastq, f"data/hpgp/r_repeat_summary_{sample_name}")
        os.unlink(temp_fastq)


def plot_r_repeat_distributions(summary_dir=None, output_file=None):
    """
    Reads r_repeat_summary files from the specified directory and creates distribution plots.
    
    This function:
    1. Scans the given directory for r_repeat_summary files
    2. Parses each file to extract r-repeat type counts for each sample
    3. Creates visualizations showing the distribution of r-repeats across samples:
       - A stacked bar chart showing the proportion of each r-repeat type per sample
       - A heatmap showing the relative abundance of each r-repeat type
       - Individual pie charts for each sample showing r-repeat distribution
       - Summary figures based only on copy numbers
    
    Parameters:
    -----------
    summary_dir : str or Path, optional
        Directory containing the r_repeat_summary files (default: "data/hpgp")
    output_file : str or Path, optional
        Path to save the output plots (default: "data/hpgp/r_repeat_distribution_plots.pdf")
        
    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the counts of each r-repeat type per sample
    """
    if summary_dir is None:
        summary_dir = Path("data/hpgp")
    else:
        summary_dir = Path(summary_dir)
    
    if output_file is None:
        output_file = summary_dir / "r_repeat_distribution_plots.pdf"
    
    # Find all summary files
    summary_files = list(summary_dir.glob("r_repeat_summary_*_repeat_summary.txt"))
    if not summary_files:
        print(f"No r_repeat summary files found in {summary_dir}")
        return None
    
    # Create a dictionary to hold the data
    data = {}
    
    # Process each summary file
    for summary_file in summary_files:
        # Extract sample name from filename
        sample_name = summary_file.stem.replace("r_repeat_summary_", "").replace("_repeat_summary", "")
        
        # Read the summary file
        r_repeat_counts = {}
        with open(summary_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    r_repeat_type, count = parts
                    r_repeat_counts[r_repeat_type] = int(count)
        
        data[sample_name] = r_repeat_counts
    
    # Convert to pandas DataFrame for easier plotting
    df = pd.DataFrame(data).fillna(0).astype(int)
    
    # Create a copy of the DataFrame for copy number analysis
    copy_number_data = {}
    
    # Get all unique r-repeat types
    all_r_repeats = df.index.tolist()
    
    # Define a consistent color palette for all r-repeat types
    color_palette = plt.cm.get_cmap('tab20', len(all_r_repeats))
    color_dict = {r_repeat: color_palette(i) for i, r_repeat in enumerate(all_r_repeats)}
    
    # Extract copy number information from each r-repeat type
    for sample_name in df.columns:
        copy_number_counts = {
            "0 copy": 0,
            "1 copy": 0,
            "2 copies": 0,
            "3 copies": 0,
            "4 copies": 0,
            "Other": 0
        }
        
        for r_repeat, count in df[sample_name].items():
            if "0copy" in r_repeat:
                copy_number_counts["0 copy"] += count
            elif "1copy" in r_repeat:
                copy_number_counts["1 copy"] += count
            elif "2copies" in r_repeat:
                copy_number_counts["2 copies"] += count
            elif "3copies" in r_repeat:
                copy_number_counts["3 copies"] += count
            elif "4copies" in r_repeat:
                copy_number_counts["4 copies"] += count
            else:
                copy_number_counts["Other"] += count
        
        copy_number_data[sample_name] = copy_number_counts
    
    # Convert copy number data to DataFrame
    copy_df = pd.DataFrame(copy_number_data).fillna(0).astype(int)
    
    # Calculate the total counts per sample
    sample_totals = df.sum()
    
    # Calculate proportions for detailed r-repeat types
    proportions_df = df.copy()
    for col in proportions_df.columns:
        if sample_totals[col] > 0:  # Avoid division by zero
            proportions_df[col] = proportions_df[col] / sample_totals[col]
    
    # Calculate proportions for copy number summary
    copy_proportions_df = copy_df.copy()
    for col in copy_proportions_df.columns:
        if sample_totals[col] > 0:  # Avoid division by zero
            copy_proportions_df[col] = copy_proportions_df[col] / sample_totals[col]
    
    # --- DETAILED R-REPEAT TYPE PLOTS ---
    
    # Create the overview plots
    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 16))
    
    # 1. Stacked bar chart of proportions
    proportions_df.transpose().plot(kind='bar', stacked=True, ax=ax1, color=[color_dict[r] for r in df.index])
    ax1.set_title('Proportion of r-repeat types by sample', fontsize=16)
    ax1.set_xlabel('Sample', fontsize=14)
    ax1.set_ylabel('Proportion', fontsize=14)
    ax1.legend(title='r-repeat type', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Annotate bars with total counts
    for i, (sample, total) in enumerate(sample_totals.items()):
        ax1.text(i, 1.01, f"n={total}", ha='center', fontsize=9)
    
    # 2. Heatmap of counts - normalize by column (sample) maximum
    norm_data = df.copy()
    for col in norm_data.columns:
        col_max = norm_data[col].max()
        if col_max > 0:  # Avoid division by zero
            norm_data[col] = norm_data[col] / col_max
    
    sns.heatmap(norm_data, annot=df, fmt="d", cmap="YlGnBu", ax=ax2, 
                cbar_kws={'label': 'Normalized Count (by sample)'})
    ax2.set_title('Number of reads per r-repeat type (normalized within each sample)', fontsize=16)
    ax2.set_xlabel('Sample', fontsize=14)
    ax2.set_ylabel('r-repeat type', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close(fig1)
    
    # --- COPY NUMBER SUMMARY PLOTS ---
    
    # Create copy number summary plots
    fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
    
    # Define consistent colors for copy numbers
    copy_colors = {
        "0 copy": "#a9a9a9",    # gray
        "1 copy": "#1f77b4",    # blue
        "2 copies": "#ff7f0e",  # orange
        "3 copies": "#2ca02c",  # green
        "4 copies": "#d62728",  # red
        "Other": "#9467bd"      # purple
    }
    
    # 1. Stacked bar chart of copy number proportions
    copy_proportions_df.transpose().plot(kind='bar', stacked=True, ax=ax1, 
                                        color=[copy_colors[c] for c in copy_df.index])
    ax1.set_title('Proportion of r-repeat copies by sample', fontsize=16)
    ax1.set_xlabel('Sample', fontsize=14)
    ax1.set_ylabel('Proportion', fontsize=14)
    ax1.legend(title='Copy number', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Annotate bars with total counts
    for i, (sample, total) in enumerate(sample_totals.items()):
        ax1.text(i, 1.01, f"n={total}", ha='center', fontsize=9)
    
    # 2. Heatmap of copy number counts
    norm_copy_data = copy_df.copy()
    for col in norm_copy_data.columns:
        col_max = norm_copy_data[col].max()
        if col_max > 0:
            norm_copy_data[col] = norm_copy_data[col] / col_max
    
    sns.heatmap(norm_copy_data, annot=copy_df, fmt="d", cmap="YlGnBu", ax=ax2,
                cbar_kws={'label': 'Normalized Count (by sample)'})
    ax2.set_title('Number of reads per copy number (normalized within each sample)', fontsize=16)
    ax2.set_xlabel('Sample', fontsize=14)
    ax2.set_ylabel('Copy number', fontsize=14)
    
    plt.tight_layout()
    copy_summary_file = Path(str(output_file).replace('.pdf', '_copy_summary.pdf'))
    plt.savefig(copy_summary_file)
    plt.close(fig3)
    
    # --- INDIVIDUAL PIE CHARTS ---
    
    # Create individual pie charts for each sample
    n_samples = len(df.columns)
    n_cols = min(3, n_samples)  # Maximum 3 plots per row
    n_rows = (n_samples + n_cols - 1) // n_cols  # Ceiling division
    
    # Detailed r-repeat type pie charts
    fig2, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    if n_samples == 1:
        axes = np.array([axes])  # Ensure axes is always iterable
    axes = axes.flatten()
    
    for i, (sample, data_col) in enumerate(df.items()):
        # Filter zero values for cleaner pie charts
        filtered_data = data_col[data_col > 0]
        
        # Create pie chart with consistent colors
        ax = axes[i]
        wedges, texts, autotexts = ax.pie(
            filtered_data,
            autopct='%1.1f%%',
            startangle=90,
            colors=[color_dict[idx] for idx in filtered_data.index],
            wedgeprops={'linewidth': 1, 'edgecolor': 'white'}
        )
        
        # Improve text contrast
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontsize(9)
        
        # Add title and legend
        ax.set_title(f"Sample: {sample}\nTotal reads: {sample_totals[sample]}", fontsize=12)
        ax.legend(filtered_data.index, title="R-repeat Types", 
                 loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    
    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
        
    plt.tight_layout()
    pie_chart_file = Path(str(output_file).replace('.pdf', '_pie_charts.pdf'))
    plt.savefig(pie_chart_file)
    plt.close(fig2)
    
    # Copy number pie charts
    fig4, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    if n_samples == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for i, (sample, data_col) in enumerate(copy_df.items()):
        # Filter zero values
        filtered_data = data_col[data_col > 0]
        
        # Create pie chart with consistent colors
        ax = axes[i]
        wedges, texts, autotexts = ax.pie(
            filtered_data,
            autopct='%1.1f%%',
            startangle=90,
            colors=[copy_colors[idx] for idx in filtered_data.index],
            wedgeprops={'linewidth': 1, 'edgecolor': 'white'}
        )
        
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontsize(9)
        
        ax.set_title(f"Sample: {sample}\nTotal reads: {sample_totals[sample]}", fontsize=12)
        ax.legend(filtered_data.index, title="Copy Numbers", 
                 loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    
    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
        
    plt.tight_layout()
    copy_pie_chart_file = Path(str(output_file).replace('.pdf', '_copy_pie_charts.pdf'))
    plt.savefig(copy_pie_chart_file)
    plt.close(fig4)
    
    print(f"Overview plots saved to {output_file}")
    print(f"Copy number summary plots saved to {copy_summary_file}")
    print(f"Individual r-repeat pie charts saved to {pie_chart_file}")
    print(f"Individual copy number pie charts saved to {copy_pie_chart_file}")
    
    return df, copy_df


def plot_r_repeat_distributions_large(summary_dir=None, output_file=None, threshold=50):
    """
    Creates a simplified plot of r-repeat copy number proportions by sample,
    with hierarchical clustering and a subtle color palette.
    
    This function:
    1. Scans the given directory for r_repeat_summary files
    2. Performs hierarchical clustering on samples based on distribution similarity
    3. Creates a visualization without individual sample names
    4. Uses a subtle color palette for better visual appeal
    
    Parameters:
    -----------
    summary_dir : str or Path, optional
        Directory containing the r_repeat_summary files (default: "data/hpgp")
    output_file : str or Path, optional
        Path to save the output plot (default: "data/hpgp/r_repeat_proportions_clustered.pdf")
    threshold : int, optional
        Minimum number of samples to use the large-scale visualization approach (default: 50)
        
    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the copy number counts per sample
    """
    if summary_dir is None:
        summary_dir = Path("data/hpgp")
    else:
        summary_dir = Path(summary_dir)
    
    if output_file is None:
        output_file = summary_dir / "r_repeat_proportions_clustered.pdf"
    
    # Find all summary files
    summary_files = list(summary_dir.glob("r_repeat_summary_*_repeat_summary.txt"))
    if not summary_files:
        print(f"No r_repeat summary files found in {summary_dir}")
        return None
        
    print(f"Processing {len(summary_files)} samples...")
    
    # Load data from summary files
    data = {}
    for summary_file in summary_files:
        # Extract sample name from filename
        sample_name = summary_file.stem.replace("r_repeat_summary_", "").replace("_repeat_summary", "")
        
        # Read the summary file
        r_repeat_counts = {}
        with open(summary_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    r_repeat_type, count = parts
                    r_repeat_counts[r_repeat_type] = int(count)
        
        data[sample_name] = r_repeat_counts
    
    # Convert to pandas DataFrame for easier analysis
    df = pd.DataFrame(data).fillna(0).astype(int)
    
    # Create copy number summary data - add 0 copy
    copy_df = pd.DataFrame(index=["0 copy", "1 copy", "2 copies", "3 copies", "4 copies", "Other"], 
                         columns=df.columns).fillna(0)
    
    # Extract copy number information
    for sample_name in df.columns:
        for r_repeat, count in df[sample_name].items():
            if "0copy" in r_repeat:
                copy_df.loc["0 copy", sample_name] += count
            elif "1copy" in r_repeat:
                copy_df.loc["1 copy", sample_name] += count
            elif "2copies" in r_repeat:
                copy_df.loc["2 copies", sample_name] += count
            elif "3copies" in r_repeat:
                copy_df.loc["3 copies", sample_name] += count
            elif "4copies" in r_repeat:
                copy_df.loc["4 copies", sample_name] += count
            else:
                copy_df.loc["Other", sample_name] += count
    
    # Calculate total reads per sample
    sample_totals = df.sum()
    
    # Calculate proportions
    copy_proportions_df = copy_df.copy()
    for col in copy_proportions_df.columns:
        if sample_totals[col] > 0:
            copy_proportions_df[col] = copy_proportions_df[col] / sample_totals[col]
    
    # ----- Perform hierarchical clustering on samples -----
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    
    # Transpose for clustering by sample
    copy_proportions_array = copy_proportions_df.fillna(0).T.values
    
    # Compute linkage matrix using Ward's method (minimizes variance within clusters)
    sample_linkage = linkage(copy_proportions_array, method='ward', optimal_ordering=True)
    
    # Determine suitable number of clusters
    n_clusters = max(3, min(20, int(np.log2(len(df.columns)))))
    sample_clusters = fcluster(sample_linkage, n_clusters, criterion='maxclust')
    
    # Create cluster assignments
    sample_to_cluster = {df.columns[i]: sample_clusters[i] for i in range(len(df.columns))}
    
    # ----- Create a dendrogram figure to show clustering hierarchy -----
    plt.figure(figsize=(12, 8))
    
    # Plot the dendrogram showing hierarchical structure without labels
    dendrogram(
        sample_linkage,
        labels=None,  # No sample labels
        leaf_rotation=0,
        color_threshold=0.7*max(sample_linkage[:,2]),  # Color threshold to highlight clusters
        above_threshold_color='gray'
    )
    
    plt.title('Hierarchical Clustering of Samples by r-repeat Copy Number Patterns', fontsize=14)
    plt.xlabel('Samples (indexed by cluster order)', fontsize=12)
    plt.ylabel('Distance', fontsize=12)
    plt.tight_layout()
    dendrogram_file = Path(str(output_file).replace('.pdf', '_dendrogram.pdf'))
    plt.savefig(dendrogram_file, dpi=300)
    plt.close()
    
    # ----- Get clustered order of samples from the dendrogram -----
    # Get the order of leaf nodes from the dendrogram and use this to order samples
    new_col_order = [df.columns[i] for i in np.argsort(sample_clusters)]
    
    # Reorder DataFrames with clustered samples
    copy_proportions_df = copy_proportions_df[new_col_order]
    copy_df = copy_df[new_col_order]
    
    # ----- Define a subtle color palette -----
    # Pastel/muted colors that are less vivid but still distinguishable
    copy_colors = {
        "0 copy": "#999999",       # gray
        "1 copy": "#8da0cb",       # muted blue
        "2 copies": "#fc8d62",     # muted orange
        "3 copies": "#66c2a5",     # muted teal
        "4 copies": "#e78ac3",     # muted magenta
        "Other": "#a6d854"         # muted green
    }
    
    # ----- Create a heatmap visualization -----
    plt.figure(figsize=(12, 6))
    
    # Create heatmap with no x-axis labels (no sample names)
    ax = sns.heatmap(
        copy_proportions_df.T,     # Transpose for samples as rows
        cmap="YlGnBu",             # Use a more subtle colormap
        cbar_kws={'label': 'Proportion of reads'},
        linewidths=0,              # No lines between cells
        xticklabels=True,          # Show x-axis labels (copy numbers)
        yticklabels=False          # Hide y-axis labels (sample names)
    )
    
    # Add cluster separation lines and annotations
    prev_cluster = sample_clusters[0]
    cluster_boundaries = []
    
    # Find cluster boundaries
    for i, col in enumerate(new_col_order):
        current_cluster = sample_to_cluster[col]
        if current_cluster != prev_cluster:
            cluster_boundaries.append(i)
            prev_cluster = current_cluster
    
    # Draw horizontal lines at cluster boundaries
    for boundary in cluster_boundaries:
        ax.axhline(y=boundary, color='black', linestyle='-', linewidth=0.5, alpha=0.5)
    
    # Add cluster annotations on the right side
    current_cluster = sample_to_cluster[new_col_order[0]]
    start_idx = 0
    
    for i, sample in enumerate(new_col_order):
        cluster = sample_to_cluster[sample]
        if cluster != current_cluster or i == len(new_col_order) - 1:
            # Add cluster label in the middle of the cluster
            if i == len(new_col_order) - 1:  # Handle the last sample
                mid_point = (start_idx + i) / 2
                cluster_size = i - start_idx + 1
            else:
                mid_point = (start_idx + i - 1) / 2
                cluster_size = i - start_idx
            
            # Add annotation on the right side
            plt.text(
                copy_proportions_df.shape[0] + 0.5,  # Position to the right of the heatmap
                mid_point,
                f"Cluster {current_cluster}\n(n={cluster_size})",
                ha='left', va='center',
                fontsize=9, 
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=0)
            )
            
            # Update for next cluster
            current_cluster = cluster
            start_idx = i
    
    plt.title('r-repeat Copy Number Proportions (Samples Grouped by Hierarchical Clustering)', fontsize=14)
    plt.xlabel('Copy Number Type', fontsize=12)
    plt.ylabel('Samples (clustered)', fontsize=12)
    plt.xticks(rotation=0)  # No rotation needed for the copy number labels
    
    # Create a secondary figure to show the copy number stacked proportions by cluster
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    
    # Calculate average proportions for each cluster
    cluster_means = {}
    for cluster_id in range(1, n_clusters + 1):
        # Get samples in this cluster
        cluster_samples = [col for col in new_col_order if sample_to_cluster[col] == cluster_id]
        # Calculate mean proportions
        if cluster_samples:
            cluster_means[f"Cluster {cluster_id}"] = copy_proportions_df[cluster_samples].mean(axis=1)
    
    # Create DataFrame from means
    cluster_df = pd.DataFrame(cluster_means)
    
    # Plot stacked bar chart with the subtle color palette
    cluster_df.plot(
        kind='bar',
        stacked=True,
        ax=ax2,
        color=[copy_colors[idx] for idx in copy_df.index],
        width=0.7,
        alpha=0.9
    )
    
    # Add informative annotations
    for i, cluster in enumerate(cluster_df.columns):
        cluster_id = int(cluster.split()[-1])
        sample_count = list(sample_to_cluster.values()).count(cluster_id)
        ax2.text(i, 1.02, f"n={sample_count}", ha='center', fontsize=9)
    
    ax2.set_title('Average r-repeat Copy Number Distribution by Cluster', fontsize=14)
    ax2.set_xlabel('Cluster', fontsize=12)
    ax2.set_ylabel('Average Proportion', fontsize=12)
    ax2.grid(axis='y', linestyle='--', alpha=0.3)
    ax2.legend(title="Copy Numbers")
    
    # Save figures
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    
    cluster_summary_file = Path(str(output_file).replace('.pdf', '_cluster_summary.pdf'))
    fig2.savefig(cluster_summary_file, dpi=300, bbox_inches="tight")
    plt.close('all')
    
    # Create a summary CSV file with cluster assignments
    summary_data = []
    for sample_name in new_col_order:
        cluster_id = sample_to_cluster[sample_name]
        total_reads = sample_totals[sample_name]
        copy_counts = {copy_type: copy_df.loc[copy_type, sample_name] for copy_type in copy_df.index}
        copy_props = {copy_type: copy_proportions_df.loc[copy_type, sample_name] for copy_type in copy_df.index}
        
        summary_data.append({
            'Sample': sample_name,
            'Cluster': cluster_id,
            'Total_Reads': total_reads,
            **{f"{copy_type}_Count": copy_counts[copy_type] for copy_type in copy_df.index},
            **{f"{copy_type}_Prop": copy_props[copy_type] for copy_type in copy_df.index}
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_csv = Path(str(output_file).replace('.pdf', '_summary.csv'))
    summary_df.to_csv(summary_csv, index=False)
    
    print(f"Hierarchical clustering dendrogram saved to {dendrogram_file}")
    print(f"Heatmap visualization saved to {output_file}")
    print(f"Cluster summary plot saved to {cluster_summary_file}")
    print(f"Summary data with clustering saved to {summary_csv}")
    print(f"Samples organized into {n_clusters} clusters")
    
    return copy_df

if __name__ == "__main__":
    # Uncomment the function you want to run
    #plot_r_repeat_distributions()
    #plot_r_repeat_distributions_large(threshold=20)  # Set threshold lower for testing
    analyze_hpgp_r_repeats()
    # make_r_repeat_refs()
    # map_r_repeats()
    # get_each_r_repeat_from_seq()