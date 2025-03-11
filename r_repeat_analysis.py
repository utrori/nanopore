from pathlib import Path
import os
import subprocess
from natsort import natsorted
import collections
import pysam
import utilities
import io


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


if __name__ == "__main__":
    analyze_hpgp_r_repeats()
    #make_r_repeat_refs()
    #map_r_repeats()
    #get_each_r_repeat_from_seq()