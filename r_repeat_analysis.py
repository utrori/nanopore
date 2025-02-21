from pathlib import Path
from natsort import natsorted
import collections
import pysam
import utilities


def make_r_repeat_refs():
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

def map_r_repeats():
    fq = "test_files/dorado_output_PSCA0047/PSCA0047_dorado.fastq"
    refs = Path('references/r_repeats').glob('*')
    summary = {} # key: read_id, value: dict with keys as r_repeat ref names and values as the alignment score
    r_repeats = []
    for ref in refs:
        r_repeat = ref.stem
        r_repeats.append(r_repeat)
        sam = utilities.minimap2_mapping(fq, ref)
        with pysam.AlignmentFile(sam, "r") as f:
            for read in f:
                if not read.is_mapped or read.query_length < 30000:
                    continue
                alignment_score = read.get_tag("AS")
                if alignment_score < 500:
                    continue
                if read.query_name not in summary:
                    summary[read.query_name] = {}
                summary[read.query_name][r_repeat] = alignment_score
    each_read_r_repeat_ret = ''
    read_id2max_r_repeat = {}
    for read_id, scores in natsorted(summary.items()):
        max_r_repeat = max(scores, key=lambda k: scores[k])
        #print(read_id, max_r_repeat)
        read_id2max_r_repeat[read_id] = max_r_repeat
        each_read_r_repeat_ret += f'{read_id}\t{max_r_repeat}\n'
    with open('data/read_id2r_repeat_type.txt', 'w') as fw:
        fw.write(each_read_r_repeat_ret)
    max_dist = collections.defaultdict(int) 
    for r_repeat in read_id2max_r_repeat.values():
        max_dist[r_repeat] += 1
    summary_ret = ''
    for r_repeat in sorted(max_dist):
        summary_ret += f'{r_repeat}\t{max_dist[r_repeat]}\n'
    with open('data/r_repeat_summary.txt', 'w') as fw:
        fw.write(summary_ret)


def get_each_r_repeat_from_seq():
    #fq = "test_files/dorado_output_PSCA0047/PSCA0047_dorado.fastq"
    fq = "test_files/201020_47.fastq"
    consensus_5 = 'TGAGACTCAGCCGGCGTCTCGCCGTGTCCCGGGTCGACCGGCGGGCCTTCTCCACCGAGCGGCGTGTAGGAGTGCCCGTCGGGACGAACCGCAACCGGAGCGTCCCCGTCTCGGTCGGCACCTCCGGGGTCGACCAGCTGCCGCCCGCGAGCTCCGGACTTAGCCGGCGCCTGCACGTGTCCCGGGTCGACCAGCAGGCGGCCGCCGGACGCTGCGGCGCACCGACGCGAGGGCGTCGATTCCCGTTCGCGCGCCCGCGACCTCCACCGGCCTCGGCCCGCGGTGGAGCTGGGACCACGCGGAACTCCCTCTCTCACATTTTTTTCAGCCCCACCGCGAGTTTGCGTCCGCGGGACTTTTAAGAGGGAGTCACTGCTGCCGTCAGCCAGTAATGCTTCCTCCTTTTTTGCTTTT'
    consensus_3 = 'TCCTTGGTGCCTTCTCGGCTC'
    consensus_5_file = 'references/r_repeat_5_consensus.fa'
    consensus_3_file = 'references/r_repeat_3_consensus.fa'
    with open(consensus_5_file, 'w') as fw:
        fw.write(f'>r_repeat_5_consensus\n{consensus_5}')
    with open(consensus_3_file, 'w') as fw:
        fw.write(f'>r_repeat_3_consensus\n{consensus_3}')
    sam = utilities.minimap2_mapping(fq, consensus_5_file)
    sam2 = utilities.bwa_mapping(fq, consensus_3_file)
    read_id2alignments = collections.defaultdict(list)
    with pysam.AlignmentFile(sam2, "r") as f:
        for read in f:
            if read.is_mapped:
                read_id2alignments[read.query_name].append((read.reference_start, read.reference_end))
    for rid, alignments in read_id2alignments.items():
        print(alignments)


if __name__ == "__main__":
    #make_r_repeat_refs()
    #map_r_repeats()
    get_each_r_repeat_from_seq()