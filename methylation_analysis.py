import numpy as np
import pickle
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


def methylation_binning(window: int, methylation_data: list[tuple[int, int]], method: str='average'):
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
    methylation_in_ref = {pos: None for pos in ref_cpg_positions}
    
    # Map read positions to reference positions for CpG sites
    for read_pos, ref_pos in alignment.get_aligned_pairs():
        if ref_pos in ref_cpg_positions and read_pos is not None:
            prob = methylation_dict.get(read_pos, -1)
            # Skip positions where CpG was present in reference but not correctly called in the read
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
    

def make_methylation_summary_in_reference_positions_with_aligned_bam(aligned_bam, called_bam, ref):
    """
    Analyze methylation data in a BAM file and summarize it by reference CpG position for each read.
    
    Parameters:
        aligned_bam: Path to BAM file with aligned reads
        called_bam: Path to BAM file with methylation calls
        ref: Path to reference sequence FASTA file
    Returns:
        Tuple containing:
        - read_methylation_dict: Dictionary with structure:
          {read_id: {
              'alignments': [
                  {
                      'reference_start': start_pos,
                      'reference_end': end_pos,
                      'query_start': q_start,
                      'query_end': q_end,
                      'methylation': {ref_pos: methylation_score, ...}
                  },
                  ...
              ]
          }}
        - global_methylation_summary: Dictionary mapping reference CpG positions to list of methylation scores
    """
    # Load reference sequence and CpG positions
    _, ref_cpg_positions = get_seq_and_find_cpgs(ref)
    
    # Initialize dictionaries to store methylation information
    global_methylation_summary = {pos: [] for pos in ref_cpg_positions}
    read_methylation_dict = {}
    
    # Iterate over reads in BAM file
    with dorado_bam_io.BamAnalyzer(aligned_bam, called_bam) as analyzer:
        reads_data = analyzer.analyze()
        
        for read_id, read_data in reads_data.items():
            # Initialize entry for this read
            read_methylation_dict[read_id] = {'alignments': []}
            
            # Access methylation data from the read
            methylation = read_data.methylation_data
            
            for alignment in read_data.alignments:
                # Get scores for this alignment
                scores = get_methylation_in_reference_positions(
                    alignment.aligned_segment, 
                    methylation, 
                    ref_cpg_positions
                )
                
                # Skip alignments with no methylation data
                if not any(score is not None for score in scores.values()):
                    continue
                
                # Create alignment entry with position info and methylation data
                alignment_entry = {
                    'reference_start': alignment.reference_start,
                    'reference_end': alignment.reference_end,
                    'query_start': alignment.query_alignment_start,
                    'query_end': alignment.query_alignment_end,
                    'methylation': {
                        pos: score for pos, score in scores.items() 
                        if score is not None
                    }
                }
                
                # Add to read's alignment list
                read_methylation_dict[read_id]['alignments'].append(alignment_entry)
                
                # Update global methylation summary
                for pos, score in scores.items():
                    if score is not None:
                        global_methylation_summary[pos].append(score)

    return read_methylation_dict, global_methylation_summary


def ref_methylation_hpgp():
    ref = config.RDNA_REF_HUMAN_COD
    for dorado_bam in Path('/media/owner/bbeedd59-668c-4d67-a65b-442d3272c3bd/hpgp/dorado_bc/').glob('*'):
        sample_name = dorado_bam.name
        print(f"Analyzing {sample_name}")
        mapped_bam = dorado_bam.parent.parent / 'dorado_aligned' / dorado_bam.name.replace('.bam', '_mapped.bam')
        ref_met_summary = make_methylation_summary_in_reference_positions_with_aligned_bam(mapped_bam, dorado_bam, ref)
        with open(f"pickles/hpgp_coding/{sample_name}_ref_met_summary.pkl", 'wb') as f:
            pickle.dump(ref_met_summary, f)


def methylation_per_coding(dorado_bam_path=None, mapped_bam_path=None, output_file=None):
    """
    Analyze methylation patterns in coding regions and create histogram.
    
    This function:
    1. Loads mapped BAM and methylation data from Dorado
    2. Analyzes methylation patterns in coding regions (4-16kb in length)
    3. Creates a histogram of average methylation scores
    
    Parameters:
    -----------
    dorado_bam_path : str or Path, optional
        Path to the BAM file with methylation data (default: "test_files/R10.4.1_PSCA0047.bam")
    mapped_bam_path : str or Path, optional
        Path to the BAM file with aligned reads (default: "test_files/R10.4.1_PSCA0047/R10.4.1_PSCA0047.bam")
    output_file : str or Path, optional
        Path to save the output histogram (default: "temp_figs/methylation_coding_histogram.png")
        
    Returns:
    --------
    list
        A list of average methylation scores for each read
    """
    # Set default paths if not provided
    if dorado_bam_path is None:
        dorado_bam_path = "test_files/R10.4.1_PSCA0047.bam"
    
    if mapped_bam_path is None:
        mapped_bam_path = "test_files/R10.4.1_PSCA0047/R10.4.1_PSCA0047.bam"
    
    if output_file is None:
        output_file = "temp_figs/methylation_coding_histogram.png"
    
    # Ensure the output directory exists
    output_dir = Path(output_file).parent
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Load and analyze methylation data
    with dorado_bam_io.BamAnalyzer(mapped_bam_path, dorado_bam_path) as analyzer:
        reads_data = analyzer.analyze()
    
    ave_met_scores = []
    for read_id, read_data in reads_data.items():
        if np.mean(read_data.get_quality_scores()) < 16 or read_data.has_strand_switch:
            continue
        for alignment in read_data.alignments:
            met_list = read_data.methylation_data
            in_met_scores = []
            if 4000 < alignment.length < 16000:
                for i in met_list:
                    if alignment.query_alignment_start < i[0] < alignment.query_alignment_end:
                        in_met_scores.append(i[1])
            if in_met_scores:
                ave_met_scores.append(np.mean(in_met_scores) / 255)
    
    # Create and save the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(ave_met_scores, bins=100)
    plt.title("Distribution of Methylation Scores in Coding Regions")
    plt.xlabel("Average Methylation Score")
    plt.ylabel("Number of Reads")
    plt.grid(True, alpha=0.3)
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    # Print summary statistics
    print(f"Analyzed {len(ave_met_scores)} reads with methylation in coding regions")
    if ave_met_scores:
        print(f"Average methylation score: {np.mean(ave_met_scores):.2f}")
        print(f"Results saved to: {output_file}")
    else:
        print("No methylation data found in coding regions")
    
    return ave_met_scores


def methylation_per_coding_hpgp():
    for dorado_bam in Path('/media/owner/bbeedd59-668c-4d67-a65b-442d3272c3bd/hpgp/dorado_bc/').glob('*'):
        mapped_bam = dorado_bam.parent.parent / 'dorado_aligned' / dorado_bam.name.replace('.bam', '_mapped.bam')
        sample_name = dorado_bam.name
        figname = f'temp_figs/hpgp_met_hists/{sample_name}_methylation_coding_histogram.png'
        methylation_per_coding(dorado_bam, mapped_bam, figname)


def analyze_per_base_methylation_by_read_quality():
    quality2met_probs = collections.defaultdict(list)
    # quality is binned by 1
    for dorado_bam in Path('/media/owner/bbeedd59-668c-4d67-a65b-442d3272c3bd/hpgp/dorado_bc/').glob('*'):
        print(dorado_bam)
        with pysam.AlignmentFile(dorado_bam, check_sq=False) as bam:
            for n, read in enumerate(bam):
                if n > 1000:
                    break
                if read.query_length < 10000:
                    continue
                ave_quality = np.mean(read.query_qualities)
                if read.modified_bases.get(('C', 0, 'm'), 0):
                    methylation_probs = [i[1] for i in read.modified_bases[('C', 0, 'm')]]
                    quality2met_probs[int(ave_quality)].extend(methylation_probs)
    for quality, met_probs in sorted(quality2met_probs.items()):
        plt.hist(met_probs, bins=100)
        plt.savefig(f'temp_figs/hpgp_quality_vs_per_base_prob_dist/{quality}_met_hist.png', dpi=300)
        plt.close()


def main():
    ref_methylation_hpgp()
    quit()
    """
    Example usage of methylation analysis functions.
    Demonstrates common workflows for methylation data analysis from nanopore reads.
    """
    print("Methylation Analysis Examples")
    print("============================\n")
    
    # Set up paths
    ref = config.RDNA_REF_HUMAN_COD
    dorado_bam = 'test_files/dorado_output_PSCA0047/calls_2025-02-05_T09-56-59.bam'
    aligned_bam = 'test_files/dorado_output_PSCA0047/PSCA0047_dorado_aligned.bam'
    output_dir = Path('example_output')
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Example 1: Per-read methylation data with reference mapping
    print("\n1. Analyzing per-read methylation with reference mapping")
    print("-------------------------------------------------")
    read_methylation_dict, global_summary = make_methylation_summary_in_reference_positions_with_aligned_bam(
        aligned_bam, dorado_bam, ref
    )
    
    # Display summary statistics
    cpg_count = len(global_summary)
    mapped_cpgs = sum(1 for scores in global_summary.values() if scores)
    print(f"- Reference contains {cpg_count} CpG positions")
    print(f"- Found methylation data for {mapped_cpgs} CpG positions ({mapped_cpgs/cpg_count:.1%})")
    
    # Display sample read data
    if read_methylation_dict:
        sample_read_id = next(iter(read_methylation_dict))
        read_data = read_methylation_dict[sample_read_id]
        print(f"- Sample read {sample_read_id}: {len(read_data['alignments'])} alignments")
        
        if read_data['alignments']:
            alignment = read_data['alignments'][0]
            cpg_sites = len(alignment['methylation'])
            avg_score = np.mean(list(alignment['methylation'].values())) if alignment['methylation'] else 0
            print(f"  - Alignment at ref:{alignment['reference_start']}-{alignment['reference_end']}")
            print(f"  - Contains {cpg_sites} CpG sites with average methylation score: {avg_score:.1f}")
    
    # Example 2: Analyze methylation in coding regions
    print("\n2. Analyzing methylation in coding regions")
    print("-------------------------------------------------")
    histogram_file = output_dir / 'coding_region_methylation.png'
    scores = methylation_per_coding(
        dorado_bam_path=dorado_bam,
        mapped_bam_path=aligned_bam,
        output_file=histogram_file
    )
    
    # Calculate and display statistics
    if scores:
        print(f"- Analyzed {len(scores)} reads with methylation in coding regions")
        print(f"- Average methylation score: {np.mean(scores):.3f}")
        print(f"- Methylation distribution: min={min(scores):.3f}, median={np.median(scores):.3f}, max={max(scores):.3f}")
        print(f"- Histogram saved to {histogram_file}")
    
    # Example 3: Methylation binning for visualization
    print("\n3. Methylation binning for visualization")
    print("-------------------------------------------------")
    # Get sample methylation data from the first read with data
    sample_methylation = None
    for read_id, read_data in read_methylation_dict.items():
        if read_data['alignments'] and any(len(a['methylation']) > 10 for a in read_data['alignments']):
            # Convert alignment's reference position-based methylation to query position-based
            # This is simplistic and just for demonstration
            methylation_pos = []
            for alignment in read_data['alignments']:
                for ref_pos, score in alignment['methylation'].items():
                    # Approximate query position - this is simplified
                    query_pos = alignment['query_start'] + (ref_pos - alignment['reference_start'])
                    methylation_pos.append((query_pos, score))
            
            if methylation_pos:
                sample_methylation = methylation_pos
                break
    
    if sample_methylation:
        print(f"- Sample read has {len(sample_methylation)} methylated positions")
        
        # Bin with 'average' method
        window_size = 200
        x_avg, y_avg = methylation_binning(window_size, sample_methylation, method='average')
        print(f"- Average method produced {len(x_avg)} bins of size {window_size}bp")
        
        # Bin with 'threshold' method
        x_thresh, y_thresh = methylation_binning(window_size, sample_methylation, method='threshold')
        print(f"- Threshold method produced {len(x_thresh)} bins of size {window_size}bp")
        
        # Create visualization
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Plot average method
        ax1.bar(x_avg, y_avg, width=window_size, alpha=0.7)
        ax1.set_title("Methylation Binning - Average Method")
        ax1.set_xlabel("Position in read (bp)")
        ax1.set_ylabel("Average methylation score")
        
        # Plot threshold method
        ax2.bar(x_thresh, y_thresh, width=window_size, alpha=0.7, color='orange')
        ax2.set_title("Methylation Binning - Threshold Method")
        ax2.set_xlabel("Position in read (bp)")
        ax2.set_ylabel("Threshold-based score")
        
        plt.tight_layout()
        binning_plot = output_dir / 'methylation_binning_example.png'
        plt.savefig(binning_plot)
        plt.close()
        print(f"- Binning visualization saved to {binning_plot}")
    
    print("\nAll examples completed. Output files saved to:", output_dir)

if __name__ == '__main__':
    main()