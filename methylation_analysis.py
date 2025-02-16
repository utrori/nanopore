import numpy as np
import config
import seq_analysis
import dorado_loader
import collections


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
    threshold = 255 * 0.7
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
                y.append(np.mean([1 if i > threshold else 0 for i in summary_data.get(n, [])])*10000)
            else:
                y.append(0)
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
        print(n, dinuc)
        if dinuc == 'CG':
            cpg_positions.append(n)
    return seq, cpg_positions


def find_true_cpgs(aligned_pair, ref_cpg_positions):

    quit()
    #print(aligned_pair[:20])
    

def main():
    ref = config.RDNA_REF_HUMAN
    ref_seq, ref_cpg_positions = get_seq_and_find_cpgs(ref)
    test_bam = 'test_files/201020_47/calls_2025-02-05_T09-56-59.bam'
    with dorado_loader.BamToSingleReadReader(test_bam) as bam:
        for single_read in bam:
            analyzer = seq_analysis.ReadAnalyzer(single_read, ref)
            analyzed_data = analyzer.minimap2_alignment()
            for aligned_pair in analyzed_data.aligned_pairs:
                find_true_cpgs(aligned_pair, ref_cpg_positions)

if __name__ == '__main__':
    main()