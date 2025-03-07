import numpy as np
import dorado_bam_io
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import pysam
import config
import seq_analysis
import dorado_loader
import collections


# Define the location and name of the log file
log_file_path = Path("logs/methylation_analysis.log")

# Create the directory if it doesn't exist
log_file_path.parent.mkdir(parents=True, exist_ok=True)

# Configure the logging module
logging.basicConfig(
    level=logging.DEBUG,  # Set the logging level (DEBUG, INFO, WARNING, etc.)
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",  # Time and log message format
    datefmt="%Y-%m-%d %H:%M:%S",  # Time format (year-month-day hour:minute:second)
    handlers=[
        logging.FileHandler(log_file_path, mode='w'),  # Log to a file. with 'w' the file is overwritten.
        logging.StreamHandler()  # Also log to the console
    ]
)

# Example of logging
logger = logging.getLogger(__name__)


def methylation_binning(window: int, methylation_data: list[(int, int)], method: str='average'):
    """
    Calculate binned CpG methylation values for plotting over a specified window.

    This function divides the genomic region into bins of fixed window size and computes a summary
    statistic of the methylation scores present within each bin. Two methods are supported:
    - "average": Computes the mean methylation score of all entries in a bin.
    - "threshold": Computes the proportion of scores exceeding 70% of the maximum score (i.e., 255 * 0.7)
        and scales this proportion by 10,000.

    Parameters:
            window (int): The size of the bin or window (in base pairs) to group methylation scores.
            methylation_data (list of tuples): A list where each element is a tuple (pos, score). 'pos' is the genomic
                    position (int) and 'score' is the methylation score (int).
            method (str, optional): The method for calculating the summary statistic for each bin. It must be either
                    "average" (default) to compute the mean of scores, or "threshold" to calculate the scaled proportion of scores
                    exceeding the defined threshold.

    Returns:
            tuple: Two lists (x, y) where:
                    x (list): Starting positions for each bin, calculated as bin_index * window.
                    y (list): The computed methylation value for each bin. This is either the average score or the scaled threshold
                            proportion depending on the chosen method.
    """
    met_threshold = 255 * 0.9
    unmet_threshold = 255 * 0.1
    summary_data = collections.defaultdict(list)
    x = []
    y = []
    for pos, score in methylation_data:
        summary_data[pos // window].append(score)
    max_pos = max(summary_data.keys())
    if method == 'average':
        for n in range(max_pos + 1):
            x.append(n * window)
            if summary_data.get(n, []):
                y.append(np.mean(summary_data.get(n, [])))
            else:
                y.append(0)
    elif method == 'threshold':
        for n in range(max_pos + 1):
            x.append(n * window)
            if summary_data.get(n, []):
                temp_summary = [i for i in summary_data.get(n, []) if i < unmet_threshold or met_threshold < i]
                if temp_summary:
                    y.append(np.mean([1 if i > met_threshold else 0 for i in temp_summary])*10000)
                else:
                    y.append(-2000)
            else:
                y.append(-2000)
    return x, y


def get_seq_and_find_cpgs(ref):
    """
    Extracts the nucleotide sequence from a FASTA file and identifies positions of 'CG' dinucleotides.

    Parameters:
        ref (str): Path to the FASTA file containing the reference sequence. The file is expected to have headers starting with '>'.

    Returns:
        tuple:
            - seq (str): The uppercase concatenated nucleotide sequence with headers removed.
            - cpg_positions (list of int): A list of 0-indexed positions indicating where the 'CG' dinucleotide occurs.

    Notes:
        - The function reads the file line by line, ignoring lines that begin with '>'.
        - It concatenates the remaining lines (after stripping any whitespace) and converts the full sequence to uppercase.
        - It iterates over the sequence to locate all instances of 'CG' and records their starting positions.
    """
    seq = ''
    with open(ref) as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq += line.strip()
    seq = seq.upper()
    cpg_positions = []
    for n, dinuc in [(i, seq[i:i+2]) for i in range(len(seq) - 1)]:
        if dinuc == 'CG':
            cpg_positions.append(n)
    return seq, cpg_positions


def find_true_cpgs(alignment: pysam.AlignedSegment, methylation, ref_cpg_positions):
    """
    Extracts the methylation scores that correspond to CpG sites that are present the reference sequence.
    """
    met_scores = []
    aligned_pairs = alignment.get_aligned_pairs()
    if not aligned_pairs:
        return None
    if alignment.is_reverse:
        rlen = alignment.query_length
        methylation = [(rlen - i[0] - 2, i[1]) for i in methylation]
    true_in_refs = [] 
    true_in_refs_in_ref = []
    for n, m in aligned_pairs:
        if m in ref_cpg_positions:
            true_in_refs.append(n)
            true_in_refs_in_ref.append(m)
    called_cpg_not_in_ref = []
    for pos, score in methylation:
        if pos in true_in_refs:
            met_scores.append(score)
        if pos not in true_in_refs and alignment.query_alignment_start <= pos <= alignment.query_alignment_end:
            called_cpg_not_in_ref.append(pos)
    if alignment.query_length > 30000:
        logger.info('{} {} {} {} {} {} {}'.format(alignment.query_name, alignment.is_forward, alignment.reference_start, len(met_scores)/len(called_cpg_not_in_ref) + 1, np.mean(met_scores), len(met_scores)/alignment.reference_length, alignment.query_length))
        logger.info(true_in_refs)
        logger.info(called_cpg_not_in_ref)
        #logger.info(true_in_refs_in_ref)
        #logger.info(methylation)
        #logger.info(met_scores)
    return met_scores
    

def get_methylation_in_reference_positions(alignment: pysam.AlignedSegment, methylation, ref_cpg_positions):
    """
    Maps methylation scores from a read to reference CpG positions.
    
    Parameters:
        alignment: Aligned segment from BAM file
        methylation: List of tuples (read_pos, score) with methylation data
        ref_cpg_positions: List of CpG positions in the reference
        
    Returns:
        Dictionary mapping reference CpG positions to methylation scores
    """
    # Return empty dict if no alignment data
    if not alignment.is_mapped or not alignment.get_aligned_pairs():
        return {}
    
    # Handle reverse strand reads - adjust positions
    if alignment.is_reverse:
        rlen = alignment.query_length
        methylation = [(rlen - pos - 2, score) for pos, score in methylation]
    
    # Create lookup dictionary for methylation data
    methylation_dict = {pos: score for pos, score in methylation}
    
    # Initialize dictionary with all reference CpG positions
    # Setting default value to None indicates position not covered by this read
    methylation_in_ref = {pos: None for pos in ref_cpg_positions}
    
    # Map read positions to reference positions for CpG sites
    for read_pos, ref_pos in alignment.get_aligned_pairs():
        if ref_pos in ref_cpg_positions and read_pos is not None:
            prob = methylation_dict.get(read_pos, -1)
            if prob != -1:
                methylation_in_ref[ref_pos] = prob
    
    # Optionally log detailed information for long reads
    #if alignment.query_length > 30000 and logger.isEnabledFor(logging.DEBUG):
        #covered_positions = sum(1 for score in methylation_in_ref.values() if score is not None)
        #logger.debug(f"Long read {alignment.query_name}: ref_start={alignment.reference_start}, "
                    #f"strand={'forward' if alignment.is_forward else 'reverse'}, "
                    #f"CpGs mapped={covered_positions}/{len(ref_cpg_positions)}, "
                    #f"length={alignment.query_length}")
    
    return methylation_in_ref


def make_methylation_summary_in_reference_positions(bam, ref) -> dict[int: list[int]]:
    """
    Analyze methylation data in a BAM file and summarize it by reference CpG position.
    
    Parameters:
        bam: Path to BAM file
        ref: Path to reference sequence FASTA file
    Returns:
        Dictionary mapping reference CpG positions to a list of methylation scores
    """
    # Load reference sequence and CpG positions
    _, ref_cpg_positions = get_seq_and_find_cpgs(ref)
    
    # Initialize dictionary to store methylation scores by reference position
    methylation_summary = {pos: [] for pos in ref_cpg_positions}
    
    # Iterate over reads in BAM file
    with dorado_loader.BamToSingleReadReader(bam) as bam:
        for n, single_read in enumerate(bam):
            # Extract methylation data and alignment information
            loadeddata = single_read.extract_read_data()
            methylation = loadeddata['methylation']
            analyzer = seq_analysis.ReadAnalyzer(single_read, ref)
            analyzed_data = analyzer.minimap2_alignment()
            
            # Process each alignment segment
            for alignment in analyzed_data.aligned_segments:
                if alignment.is_mapped:
                    # Get methylation scores for CpG sites in reference
                    scores = get_methylation_in_reference_positions(alignment, methylation, ref_cpg_positions)
                    # Update summary dictionary with scores
                    for pos, score in scores.items():
                        if score is not None:
                            methylation_summary[pos].append(score)
    
    return methylation_summary
    

def make_methylation_summary_in_reference_positions_with_aligned_bam(aligned_bam, called_bam, ref) -> dict[int: list[int]]:
    """
    Analyze methylation data in a BAM file and summarize it by reference CpG position.
    
    Parameters:
        aligned_bam: Path to BAM file with aligned reads
        called_bam: Path to BAM file with methylation calls
        ref: Path to reference sequence FASTA file
    Returns:
        Dictionary mapping reference CpG positions to a list of methylation scores
    """
    # Load reference sequence and CpG positions
    _, ref_cpg_positions = get_seq_and_find_cpgs(ref)
    
    # Initialize dictionary to store methylation scores by reference position
    methylation_summary = {pos: [] for pos in ref_cpg_positions}
    
    # Iterate over reads in BAM file
    with dorado_bam_io.BamAnalyzer(aligned_bam, called_bam) as analyzer:
        reads_data = analyzer.analyze()
        n = 0
        for read_id, read_data in reads_data.items():
            n += 1
            if n > 100:
                break
            for alignment in read_data.alignments:
                methylation = alignment.methylation_data
                scores = get_methylation_in_reference_positions(alignment.aligned_segment, methylation, ref_cpg_positions)
                if scores:
                    for pos, score in scores.items():
                        if score is not None:
                            methylation_summary[pos].append(score)

    return methylation_summary

def main():
    ref = config.RDNA_REF_HUMAN_COD
    ref_seq, ref_cpg_positions = get_seq_and_find_cpgs(ref)
    test_bam = 'test_files/dorado_output_PSCA0047/calls_2025-02-05_T09-56-59.bam'
    aligned_bam = 'test_files/dorado_output_PSCA0047/PSCA0047_dorado_aligned.bam'
    methylaion_summary = make_methylation_summary_in_reference_positions_with_aligned_bam(aligned_bam, test_bam, ref)
    #methylaion_summary = make_methylation_summary_in_reference_positions(test_bam, ref)
    print(methylaion_summary[12753], methylaion_summary[241])
    quit()
    all_met_scores = []
    true_met_scores = []
    num_reads_per_sample = 1000
    phred_scores = []
    with dorado_loader.BamToSingleReadReader(test_bam) as bam:
        logger.info(f'Analyzing {test_bam}')
        for n, single_read in enumerate(bam):
            loadeddata = single_read.extract_read_data()
            methylation = loadeddata['methylation']
            average_phred = np.mean(loadeddata['quality_scores'])
            phred_scores.append(average_phred)
            if average_phred < 15:
                continue
            analyzer = seq_analysis.ReadAnalyzer(single_read, ref)
            analyzed_data = analyzer.minimap2_alignment()
            for alignment in analyzed_data.aligned_segments:
                if alignment.is_mapped:
                    scores = get_methylation_in_reference_positions(alignment, methylation, ref_cpg_positions)
                    print(scores)
                    #met_scores = [i[1] for i in methylation]
                    #if scores:
                        #all_met_scores.extend(met_scores)
                        #true_met_scores.extend(scores)
            if n > num_reads_per_sample:
                break
    quit()
    mid_ratio_all = len([i for i in all_met_scores if 30 < i < 220])/len(all_met_scores)
    mid_ratio_true = len([i for i in true_met_scores if 30 < i < 220])/len(true_met_scores)
    print(mid_ratio_all, mid_ratio_true)
    quit()
    plt.hist(all_met_scores, bins=100)
    plt.savefig('./temp_figs/all_met_scores.png')
    plt.close()
    plt.hist(true_met_scores, bins=100)
    plt.savefig('./temp_figs/true_met_scores.png')


if __name__ == '__main__':
    main()