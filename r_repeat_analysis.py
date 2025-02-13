from pathlib import Path
import pysam
import subprocess
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
    with open("references/r_repeats/r_repeat_4copies_1233.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat1_seq}{repeat2_seq}{repeat3_seq}{repeat3_seq}{seq_3}\n")
    with open("references/r_repeats/r_repeat_3copies.fa", 'w') as f:
        f.write(f">r_repeat_3copies\n{seq_5}{repeat1_seq}{repeat2_seq}{repeat3_seq}{seq_3}\n")
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
    fq = "test_files/201020_47/201020_47.fastq"
    r_repeats = ["r_repeat_4copies_1123", "r_repeat_4copies_1233", "r_repeat_3copies", "r_repeat_2copies_12", "r_repeat_2copies_13", "r_repeat_2copies_23", "r_repeat_1copy_1", "r_repeat_1copy_2", "r_repeat_1copy_3"]
    summary = {} # key: read_id, value: dict with keys as r_repeat ref names and values as the alignment score
    for r_repeat in r_repeats:
        sam = utilities.minimap2_mapping(fq, f"references/r_repeats/{r_repeat}.fa")
        with pysam.AlignmentFile(sam, "r") as f:
            for read in f:
                if not read.is_mapped or read.query_length < 30000:
                    continue
                alignment_score = read.get_tag("AS")
                if alignment_score < 300:
                    continue
                if read.query_name not in summary:
                    summary[read.query_name] = {}
                summary[read.query_name][r_repeat] = alignment_score
    for read_id, scores in summary.items():
        print(read_id)
        for r_repeat, score in scores.items():
            print(f"{r_repeat}: {score}")
        print("-" * 20)
    



if __name__ == "__main__":
    make_r_repeat_refs()
    map_r_repeats()