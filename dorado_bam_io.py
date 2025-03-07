import pysam
import matplotlib.pyplot as plt
import pathlib as Path
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional

@dataclass
class AlignmentInfo:
    reference_name: str
    reference_start: int
    reference_end: int
    query_alignment_start: int
    query_alignment_end: int
    length: int
    cigar: str
    mapping_quality: int
    flag: int
    methylation_data: List[Tuple[int, int]] = field(default_factory=list)
    aligned_segment: pysam.AlignedSegment = None

@dataclass
class ReadData:
    read_id: str
    alignments: List[AlignmentInfo] = field(default_factory=list)

    def add_alignment(self, alignment: AlignmentInfo):
        self.alignments.append(alignment)

class BamAnalyzer:
    # Mapped bam needs to be cleaned with samtools view -F 4
    def __init__(self, mapped_bam_path: str, dorado_bam_path: str):
        self.mapped_bam_path = mapped_bam_path
        self.dorado_bam_path = dorado_bam_path
        self.mapped_bam = None
        self.dorado_bam = None

    def __enter__(self):
        self.dorado_bam = pysam.AlignmentFile(self.dorado_bam_path, "rb", check_sq=False)
        self.mapped_bam = pysam.AlignmentFile(self.mapped_bam_path, "rb", check_sq=False)
        self.methylation_data_by_read_id, self.quality_data_by_read_id = self._load_dorado_data()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.mapped_bam.close()
        self.dorado_bam.close()

    def _load_dorado_data(self) -> Dict[str, List[Tuple[int, int]]]:
        """Loads CpG methylation data from the Dorado BAM into a dictionary."""
        methylation_data_by_read_id = {}
        quality_by_read_id = {}
        
        for read in self.dorado_bam:
            read_id = read.query_name
            if read.has_tag('MM') or read.has_tag('Ml'):
                modified_bases = read.modified_bases
                if ('C', 0, 'm') in modified_bases:
                    methylation_data_by_read_id[read_id] = read.modified_bases[('C', 0, 'm')]
            quality = read.query_qualities
            quality_by_read_id[read_id] = quality

        return methylation_data_by_read_id, quality_by_read_id

    def analyze(self) -> Dict[str, ReadData]:
        """Analyzes multiply-aligned reads, extracting methylation data and organizing by read ID."""
        reads_data = {}

        for read in self.mapped_bam:

            read_id = read.query_name

            if read_id not in reads_data:
                reads_data[read_id] = ReadData(read_id)
            
            if read.is_supplementary:
                hard_clipped = read.cigartuples[-1][1]
            else:
                hard_clipped = 0

            # Calculate the query alignment start and end based on the read's orientation
            query_len = read.query_length + hard_clipped
            if read.is_reverse:
                query_alignment_start = query_len - read.query_alignment_end
                query_alignment_end = query_len - read.query_alignment_start
            else:
                query_alignment_start = read.query_alignment_start
                query_alignment_end = read.query_alignment_end

            alignment_info = AlignmentInfo(
                reference_name=self.mapped_bam.get_reference_name(read.reference_id),
                reference_start=read.reference_start,
                reference_end=read.reference_end,
                query_alignment_start=query_alignment_start,
                query_alignment_end=query_alignment_end,
                length = read.reference_end - read.reference_start,
                cigar=read.cigarstring,
                mapping_quality=read.mapping_quality,
                flag=read.flag,
                methylation_data=self.methylation_data_by_read_id.get(read_id, []),
                aligned_segment=read
            )
            reads_data[read_id].add_alignment(alignment_info)

        return reads_data

def methylation_per_coding():
    dorado_bam_path = "test_files/R10.4.1_PSCA0047.bam"
    mapped_bam_path = "test_files/R10.4.1_PSCA0047/R10.4.1_PSCA0047.bam"
    #dorado_bam_path = "test_files/220705_Y9993fN/calls_2025-02-05_T07-37-34.bam"
    #mapped_bam_path = "test_files/220705_Y9993fN_coding_mapped.bam"
    with BamAnalyzer(mapped_bam_path, dorado_bam_path) as analyzer:
        reads_data = analyzer.analyze()
    ave_met_scores = []
    for read_id, read_data in reads_data.items():
        for alignment in read_data.alignments:
            met_list = alignment.methylation_data
            in_met_scores = []
            if 4000 < alignment.length < 16000:
                for i in met_list:
                    if alignment.query_alignment_start < i[0] < alignment.query_alignment_end:
                        in_met_scores.append(i[1])
            if in_met_scores:
                ave_met_scores.append(np.mean(in_met_scores))
    plt.hist(ave_met_scores, bins=100)
    plt.show()


def main():
    # Example of how to use the BamAnalyzer class
    methylation_per_coding()
    quit()
    dorado_bam_path = "test_files/201020_47/calls_2025-02-05_T09-56-59.bam"
    mapped_bam_path = "test_files/201020_47_coding_mapped.bam"

    with BamAnalyzer(mapped_bam_path, dorado_bam_path) as analyzer:
      reads_data = analyzer.analyze()

    # Example of how to access and process the results
    for read_id, read_data in reads_data.items():
        print(f"Read ID: {read_id}")
        for alignment in read_data.alignments:
            print(f"  Reference: {alignment.reference_name}:{alignment.reference_start}-{alignment.reference_end}")
            print(f"  Query Alignment: {alignment.query_alignment_start}-{alignment.query_alignment_end}")
            print(f"  CIGAR: {alignment.cigar}")
            print(f"  MAPQ: {alignment.mapping_quality}")
            print(f"  Flag: {alignment.flag}")
            if alignment.methylation_data:
                for meth_tag, ml_tag in alignment.methylation_data:
                    print(f"    Modified Bases (MM): {meth_tag}")
                    print(f"    Modification Probabilities (ML): {ml_tag}")
            else:
                print("    No methylation data found.")
        print("-" * 20)

if __name__ == "__main__":
    main()