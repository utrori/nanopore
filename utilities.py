import subprocess
import os
import tempfile
import pysam
from pathlib import Path
import numpy as np
import matplotlib.collections as mc
import matplotlib.pyplot as plt
import io

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


def make_temp_full_fastq(header, seq, quality):
    """Creates a temporary FASTQ file for minimap2 alignment.

    Args:
        header (str): The FASTQ header (without the leading '@').
        read (str): The read sequence.
        quality (str): The quality string.
    """

    with tempfile.NamedTemporaryFile(mode='w+t', suffix=".fastq", delete=False, prefix="temp_fastq_") as temp_file:
        temp_file_path = Path(temp_file.name)
        fw = temp_file
        temp_file.write(f"@{header}\n{seq}\n+\n{quality}\n")
        temp_file.close()
    return temp_file_path

def minimap2_mapping(in_fastq, ref_path, out_bam=None, only_mapped=False):
    """Performs minimap2 mapping of reads in a FASTQ file to a reference genome.

    Args:
        in_fastq (Path or str): Path to the input FASTQ file.
        ref_path (Path or str): Path to the reference genome FASTA file.
        out_bam (Path or str, optional): Path to output BAM file. If provided, mapping results will be saved here.
        only_mapped (bool, optional): If True, only mapped reads will be included in the output BAM. Default is False.

    Returns:
        Path or None: Path to the temporary SAM file if out_bam is None, otherwise None.
    """
    # Convert paths to strings for subprocess if they are Path objects
    in_fastq_str = str(in_fastq)
    ref_path_str = str(ref_path)
    
    # Basic minimap2 command
    minimap_cmd = ["minimap2", "-ax", "map-ont", "-t", "16", "-Y", ref_path_str, in_fastq_str]
    
    if out_bam:
        # Direct pipe to samtools for output to BAM
        minimap_process = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE)
        
        # Construct samtools command based on whether we want only mapped reads
        samtools_cmd = ["samtools", "view", "-b"]
        if only_mapped:
            samtools_cmd.extend(["-F", "4"])  # Filter out unmapped reads
            
        samtools_cmd.extend(["-o", str(out_bam), "-"])
        
        # Execute the samtools command with input from minimap2
        subprocess.run(samtools_cmd, stdin=minimap_process.stdout, check=True)
        minimap_process.stdout.close()
        minimap_process.wait()
        
        return None
    else:
        # If no output BAM is specified, save to temporary SAM file and return its path
        result = subprocess.run(minimap_cmd, capture_output=True, text=True, check=True)
        sam_output = result.stdout
        
        with tempfile.NamedTemporaryFile(mode='w+t', suffix=".sam", delete=False, prefix="temp_sam_") as temp_file:
            temp_file.write(sam_output)
            return Path(temp_file.name)

def minimap_mapping_and_sam_analysis(header, read, quality, ref_path):
    """Maps a read using minimap2 and parses the SAM output.

    Args:
        header (str): The FASTQ header.
        read (str): The read sequence.
        quality (str): The quality string.
        ref_path (Path or str): Path to the reference genome FASTA file.

    Returns:
        mapped_regions (list): A list of tuples containing the mapped regions.
        cigar_converted (list): A list of CIGAR strings for each alignment.
    """
    
    temp_fastq_path = make_temp_full_fastq(header, read, quality)
    temp_sam_path = minimap2_mapping(temp_fastq_path, ref_path)

    mapped_regions = []
    coneverted_cigars = []
    alignments = []
    with pysam.AlignmentFile(temp_sam_path, 'r') as samfile:
        for alignment in samfile:
            if alignment.is_mapped:
                mapped_regions.append((alignment.reference_start, alignment.reference_end))
                coneverted_cigars.append(alignment.cigarstring)
                alignments.append(alignment)
    # Clean up temporary files
    temp_fastq_path.unlink()
    temp_sam_path.unlink()

    return mapped_regions, coneverted_cigars, alignments

def make_temp_fastq(header, read, quality, split_length):
    """Creates a temporary FASTQ file by splitting a read into smaller chunks.

    Args:
        header (str): The FASTQ header (without the leading '@').
        read (str): The read sequence.
        quality (str): The quality string.
        split_length (int): The length to split the read into.
    Returns:
        temp_file_path (Path): Path to the temporary FASTQ file.
    """
    split_reads = split_sequence(read, split_length)
    split_qualities = split_sequence(quality, split_length)
    
    with tempfile.NamedTemporaryFile(mode='w+t', suffix=".fastq", delete=False, prefix="temp_fastq_") as temp_file:
        temp_file_path = Path(temp_file.name)
        fw = temp_file
        for i, (split_read, split_quality) in enumerate(zip(split_reads, split_qualities)):
            fw.write(f"@{header}_{i+1}\n{split_read}\n+\n{split_quality}\n")
        temp_file.close()
    return temp_file_path

def bwa_mapping(in_fastq, ref_path, out_bam=None, only_mapped=False, multi=False):
    """Performs BWA mapping of reads in a FASTQ file to a reference genome.

    Args:
        in_fastq (Path or str): Path to the input FASTQ file.
        ref_path (Path or str): Path to the reference genome FASTA file.
        out_bam (Path or str, optional): Path to output BAM file. If provided, mapping results will be saved here.
        only_mapped (bool, optional): If True, only mapped reads will be included in the output BAM. Default is False.
        multi (bool, optional): If True, use the -a option for multiple hits (suitable for short reads)
                      otherwise, use the -M option for marking shorter split hits as secondary
    
    Returns:
        Path or None: Path to the temporary SAM file if out_bam is None, otherwise None.
    """
    # Convert paths to strings for subprocess if they are Path objects
    ref_path_str = str(ref_path)
    in_fastq_str = str(in_fastq)

    # Build the BWA command based on multi parameter
    if multi:
        bwa_cmd = ["bwa", "mem", "-Ma", "-x", "ont2d", "-t", "10", "-Y", ref_path_str, in_fastq_str]
    else:
        bwa_cmd = ["bwa", "mem", "-M", "-x", "ont2d", "-t", "10", "-Y", ref_path_str, in_fastq_str]
    
    if out_bam:
        # Direct pipe to samtools for output to BAM
        bwa_process = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
        
        # Construct samtools command based on whether we want only mapped reads
        samtools_cmd = ["samtools", "view", "-b"]
        if only_mapped:
            samtools_cmd.extend(["-F", "4"])  # Filter out unmapped reads
            
        samtools_cmd.extend(["-o", str(out_bam), "-"])
        
        # Execute the samtools command with input from bwa
        subprocess.run(samtools_cmd, stdin=bwa_process.stdout, check=True)
        bwa_process.stdout.close()
        bwa_process.wait()
        
        return None
    else:
        # If no output BAM is specified, save to temporary SAM file and return its path
        result = subprocess.run(bwa_cmd, capture_output=True, text=True, check=True)
        sam_output = result.stdout
        
        with tempfile.NamedTemporaryFile(mode='w+t', suffix=".sam", delete=False, prefix="temp_sam_") as temp_file:
            temp_file.write(sam_output)
            return Path(temp_file.name)

def split_mapping_and_sam_analysis(split_length: int, header: str, seq: str, quality: str, ref_path: str, out_bam=None):
    """Splits a read, maps it using BWA, and parses the SAM output.

    Args:
        split_length (int): Length to split the read into.
        header (str): Header of the read.
        seq (str): Read sequence.
        quality (str): Quality string.
        ref_path (Path or str): Path to the reference genome FASTA file.
        out_bam (Path or str, optional): If provided, save the mapping results to this BAM file.

    Returns:
        np.ndarray: A NumPy array containing the parsed SAM records.
    """
    
    temp_fastq_path = make_temp_fastq(header, seq, quality, split_length)
    
    if out_bam:
        # If out_bam is provided, save the mapping results directly
        bwa_mapping(temp_fastq_path, ref_path, out_bam=out_bam)
        # We still need to parse the SAM for analysis, so create a temporary SAM file
        temp_sam_path = bwa_mapping(temp_fastq_path, ref_path)
    else:
        temp_sam_path = bwa_mapping(temp_fastq_path, ref_path)

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
            # the structure of sam_info is [flag, direction, position, CIGAR, AS score, read_length]

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

def get_read_structure(read_id, split_length, samdata, species, offset=0, savename=None, title=None):
    """Plots the structure of a read based on its split-mapped SAM data.

    Args:
        read_id (str): The read ID.
        split_length (int): The length used to split the read during mapping.
        samdata (np.ndarray): The NumPy array containing the parsed SAM records.
        species (str): The species of the reference genome ('mouse', 'human', or 'fly').
        offset (int): Offset for the rDNA reference coordinates.
        savename (Path or str, optional): If provided, save the plot to this file.
        title (str, optional): Title for the plot.

    Returns:
        matplotlib.collections.LineCollection: The LineCollection object representing the plot, 
                                                or None if savename is provided.
    """
    if species == 'mouse':
        rDNA_size = 38542
        CR_size = 13403
    elif species == 'human':
        rDNA_size = 44838
        CR_size = 13314
    elif species == 'fly':
        rDNA_size = 13188
        CR_size = 13314

    plot_data = []
    for flag, _, pos, _, _, _ in samdata:
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


def bam_concatenate(bam_files: list[str], output_bam: str):
    """Concatenates multiple BAM files into a single BAM file."""
    bam_strs = [str(bam) for bam in bam_files] # Convert Path objects to strings

    command = f"samtools cat -@ 10 -o {output_bam} {' '.join(bam_strs)}"
    print(command)
    run_command(command, verbose=True)


def dorado_basecalling(pod5_file: str, output_bam: str):
    command = f'~/Softwares/dorado-0.9.1-linux-x64/bin/dorado basecaller --modified-bases-models ~/Softwares/dorado-0.9.1-linux-x64/models/dna_r9.4.1_e8_hac@v3.3_5mCG@v0.1 /home/owner/Softwares/dorado-0.9.1-linux-x64/models/dna_r9.4.1_e8_hac@v3.3 {pod5_file} > {output_bam}'
    run_command(command, verbose=True)


def dorado_bc_hpgp():
    hgpg_dir = Path("/media/owner/bbeedd59-668c-4d67-a65b-442d3272c3bd/hpgp")
    for f5 in hgpg_dir.glob('**/batch0.fast5'):
        pod5_file = f5.with_suffix('.pod5')
        command = f'pod5 convert fast5 {f5} --output {pod5_file} --threads 15'
        #run_command(command, verbose=True)
        out_bam_basedir = Path("/media/owner/bbeedd59-668c-4d67-a65b-442d3272c3bd/hpgp/dorado_bc/")
        sample = f5.parts[-2].split('_rDNA')[0]
        out_bam = out_bam_basedir / (sample + '.bam')
        #dorado_basecalling(pod5_file, out_bam)
        dorado_mapped_bam = out_bam_basedir.parent / 'dorado_aligned' / (sample + '_mapped.bam')
        temp_fastq = out_bam_basedir.parent / 'dorado_aligned' / (sample + '.fastq')
        command = f'samtools fastq -@ 10 {out_bam} > {temp_fastq}'
        #run_command(command, verbose=True)
        #minimap2_mapping(temp_fastq, "references/human_rDNA_only_coding.fa", out_bam=dorado_mapped_bam, only_mapped=True)
        os.unlink(temp_fastq)

if __name__ == "__main__":
    dorado_bc_hpgp()
    # Example usage
    sam = minimap2_mapping("test_files/HG02723_1/HG02723_1.fastq", "references/human_rDNA_only_coding.fa")
    with open(sam) as f:
        content = f.read()
    with open('test_files/HG02723_1/HG02723_1.sam', 'w') as f:
        f.write(content)