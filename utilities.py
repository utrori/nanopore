import subprocess
import pysam
import random
import os
from pathlib import Path
import numpy as np
import matplotlib.collections as mc
import matplotlib.pyplot as plt

# Use Path objects for file paths
TEMP_DIR = Path("temp_files")

# Create the temporary directory if it doesn't exist
TEMP_DIR.mkdir(parents=True, exist_ok=True)

def run_command(command, verbose=False):
    """Runs a shell command, optionally printing the command.

    Args:
        command (str): The command string to execute.
        verbose (bool): If True, print the command being executed.
    """
    if verbose:
        print(f"Running command: {command}")
    # Consider using capture_output=True to handle output more effectively.
    subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=False)

def easy_flag(flag, base):
    """Checks if a specific bit (base) is set in the SAM flag.

    Args:
        flag (int): The SAM flag value.
        base (int): The bit position to check (e.g., 1, 2, 4, 8, 16...).

    Returns:
        int: 1 if the bit is set, 0 otherwise.
    """
    return (flag // base) % 2

def circular_slice(sequence, start, end):
    """Extracts a slice from a circular sequence (e.g., a circular genome).

    Args:
        sequence (str): The sequence string.
        start (int): The starting index of the slice (can be negative).
        end (int): The ending index of the slice (can be greater than the length of the sequence).

    Returns:
        str: The sliced sequence.

    Raises:
        ValueError: If start is greater than or equal to end.
    """
    if start >= end:
        raise ValueError("Start index must be less than end index.")

    seq_len = len(sequence)
    if start < 0:
        start = seq_len + start  # Adjust for negative indexing
    if end > seq_len:
        return sequence[start:] + sequence[:end % seq_len]
    return sequence[start:end]

def split_sequence(sequence, split_length):
    """Splits a sequence into chunks of a specified length.

    Args:
        sequence (str): The sequence string.
        split_length (int): The desired length of each chunk.

    Returns:
        list: A list of sequence chunks.
    """
    return [sequence[i:i + split_length] for i in range(0, len(sequence), split_length)]

def make_temp_fastq(header, read, quality, identifier, split_length,
                    temp_file_prefix="temp_fastq"):
    """Creates a temporary FASTQ file by splitting a read into smaller chunks.

    Args:
        header (str): The FASTQ header (without the leading '@').
        read (str): The read sequence.
        quality (str): The quality string.
        identifier (str, float): A unique identifier for the read (can be a number).
        split_length (int): The length to split the read into.
        temp_file_prefix (str): Prefix for the temporary file name.
    """
    split_reads = split_sequence(read, split_length)
    split_qualities = split_sequence(quality, split_length)
    
    temp_file_path = TEMP_DIR / f"{temp_file_prefix}_{identifier}.fastq"

    with open(temp_file_path, 'w') as fw:
        for i, (split_read, split_quality) in enumerate(zip(split_reads, split_qualities)):
            fw.write(f"@{header}_{i+1}\n{split_read}\n+\n{split_quality}\n")

def bwa_mapping(ref_path, in_fastq, out_sam, multi=False):
    """Performs BWA mapping of reads in a FASTQ file to a reference genome.

    Args:
        ref_path (Path or str): Path to the reference genome FASTA file.
        in_fastq (Path or str): Path to the input FASTQ file.
        out_sam (Path or str): Path to the output SAM file.
        multi (bool): If True, use the -a option for multiple hits (suitable for short reads)
                      otherwise, use the -M option for marking shorter split hits as secondary
    """
    # Convert paths to strings for subprocess if they are Path objects
    ref_path_str = str(ref_path)
    in_fastq_str = str(in_fastq)
    out_sam_str = str(out_sam)

    if multi:
        command = ["bwa", "mem", "-Ma", "-x", "ont2d", "-t", "10", ref_path_str, in_fastq_str]
    else:
        command = ["bwa", "mem", "-M", "-x", "ont2d", "-t", "10", ref_path_str, in_fastq_str]
    
    with open(out_sam_str, "w") as outfile:
        subprocess.run(command, stdout=outfile, stderr=subprocess.DEVNULL, check=True)

def split_mapping_and_sam_analysis(split_length, header, read, quality, ref_path):
    """Splits a read, maps it using BWA, and parses the SAM output.

    Args:
        split_length (int): Length to split the read into.
        header (str): Header of the read.
        read (str): Read sequence.
        quality (str): Quality string.
        ref_path (Path or str): Path to the reference genome FASTA file.

    Returns:
        np.ndarray: A NumPy array containing the parsed SAM records.
    """
    identifier = random.random()
    
    # Use Path objects for file operations
    temp_fastq_path = TEMP_DIR / f"temp_fastq_{identifier}.fastq"
    temp_sam_path = TEMP_DIR / f"single_split_mapped_{identifier}.sam"
    
    make_temp_fastq(header, read, quality, identifier, split_length, temp_file_prefix="temp_fastq")
    bwa_mapping(ref_path, temp_fastq_path, temp_sam_path)

    sam_info = []
    with pysam.AlignmentFile(temp_sam_path, 'r') as samfile:
        for alignment in samfile:
            if alignment.is_secondary:
                continue  # Skip secondary alignments

            sam_info.append([
                alignment.flag,
                '-' if alignment.is_reverse else '+' if not alignment.is_unmapped else '0',
                alignment.reference_start,
                alignment.cigarstring,
                alignment.get_tag('AS'),
                alignment.query_length
            ])

    # Clean up temporary files
    temp_fastq_path.unlink()
    temp_sam_path.unlink()

    return np.array(sam_info)

def make_shifted_ref(ref_filepath, new_ref_filepath, shift):
    """Creates a new reference FASTA file with a shifted sequence.

    Args:
        ref_filepath (Path or str): Path to the original reference FASTA file.
        new_ref_filepath (Path or str): Path to the new shifted reference FASTA file.
        shift (int): The number of bases to shift the sequence.
    """
    
    # Convert paths to Path objects
    ref_filepath = Path(ref_filepath)
    new_ref_filepath = Path(new_ref_filepath)

    with open(ref_filepath, "r") as infile, open(new_ref_filepath, "w") as outfile:
        header = infile.readline()
        outfile.write(header)
        
        sequence = "".join(line.strip() for line in infile)
        shifted_sequence = circular_slice(sequence, shift, len(sequence) + shift)
        
        for i in range(0, len(shifted_sequence), 70):
            outfile.write(shifted_sequence[i:i + 70] + "\n")

def plot_read_structure(header, split_length, samdata, mouse=False, offset=0, savename=None, title=None):
    """Plots the structure of a read based on its split-mapped SAM data.

    Args:
        header (str): The read header.
        split_length (int): The length used to split the read during mapping.
        samdata (np.ndarray): The NumPy array containing the parsed SAM records.
        mouse (bool): If True, use mouse reference genome parameters.
        offset (int): Offset for the rDNA reference coordinates.
        savename (Path or str, optional): If provided, save the plot to this file.
        title (str, optional): Title for the plot.

    Returns:
        matplotlib.collections.LineCollection: The LineCollection object representing the plot, 
                                                or None if savename is provided.
    """
    rDNA_size = 38542 if mouse else 44838
    CR_size = 13403 if mouse else 13314

    plot_data = []
    for flag, direction, pos, _, _, _ in samdata:
        flag = int(flag)
        pos = int(pos)
        if easy_flag(flag, 4) == 1:  # Unmapped
            plot_data.append((-10000, "*"))
        else:
            adjusted_pos = (pos + offset) % rDNA_size
            plot_data.append((adjusted_pos, "+" if easy_flag(flag, 16) != 1 else "-"))

    read_num = len(plot_data)
    x_coords = np.linspace(0, split_length * (read_num - 1), read_num)
    
    vertical_lines = []
    for n, (flag, _, pos, cigar, _, read_len) in enumerate(samdata):
        flag = int(flag)
        pos = int(pos)
        read_len = int(read_len)
        if easy_flag(flag, 4) == 1:  # Unmapped
            line_segment = [(x_coords[n], -10000), (x_coords[n] + split_length, -10000)]
        elif easy_flag(flag, 16) != 1:  # Forward
            line_segment = [(x_coords[n], pos), (x_coords[n] + read_len, pos + read_len)]
        else:  # Reverse
            line_segment = [(x_coords[n], pos + read_len), (x_coords[n] + read_len, pos)]
        vertical_lines.append(line_segment)

    # Add lines for 0, CR_size, and rDNA_size
    vertical_lines.extend([
        [[0, 0], [split_length * read_num, 0]],
        [[0, CR_size], [split_length * read_num, CR_size]],
        [[0, rDNA_size], [split_length * read_num, rDNA_size]]
    ])

    # Line widths and styles
    line_widths = [1] * read_num + [0.5] * 3
    line_styles = ["solid"] * read_num + ["dashed"] * 3

    # Colors based on mapping position
    line_colors = []
    for coord, _ in plot_data:
        if coord == -10000:
            line_colors.append("black")
        else:
            line_colors.append("black")  # You can use a colormap if needed
    line_colors.extend(["black"] * 3)
    
    line_collection = mc.LineCollection(
        vertical_lines, linewidths=line_widths, linestyles=line_styles, colors=line_colors
    )

    if savename:
        # Convert savename to Path object if it's a string
        savename = Path(savename)
        
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.2)
        ax.add_collection(line_collection)
        ax.autoscale()
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
        ax.set_yticklabels(("unmapped", 0, 10000, 20000, 30000, 40000))
        ax.set_ylim((-12000, 46000))
        if title:
            ax.set_title(title)
        plt.savefig(savename, dpi=300)
        plt.close(fig)
        return None
    else:
        return line_collection

def make_rDNA_fastq(ref_filepath="rDNA_index/humRibosomal2.fa", output_fastq="rDNA_split_fastq", split_length=100):
    """Generates a FASTQ file containing split rDNA sequences.

    Args:
        ref_filepath (str): Path to the reference FASTA file (default: "rDNA_index/humRibosomal2.fa").
        output_fastq (str): Path to the output FASTQ file (default: "rDNA_split_fastq").
        split_length (int): Length to split the rDNA sequence into (default: 100).
    """
    
    # Convert paths to Path objects
    ref_filepath = Path(ref_filepath)
    output_fastq = Path(output_fastq)
    
    sequence = ""
    with open(ref_filepath, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()

    with open(output_fastq, "w") as fw:
        for n in range(0, len(sequence), split_length):
            split_seq = sequence[n:n + split_length]
            quality = "J" * len(split_seq)  # Dummy quality scores
            fw.write(f"@rDNA_{n // split_length}\n{split_seq}\n+\n{quality}\n")