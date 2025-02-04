import pysam
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
    cigar: str
    mapping_quality: int
    flag: int
    methylation_data: List[Tuple[str, list]] = field(default_factory=list)

@dataclass
class ReadData:
    read_id: str
    alignments: List[AlignmentInfo] = field(default_factory=list)

    def add_alignment(self, alignment: AlignmentInfo):
        self.alignments.append(alignment)

class BamAnalyzer:
    def __init__(self, mapped_bam_path: str, dorado_bam_path: str):
        self.mapped_bam_path = mapped_bam_path
        self.dorado_bam_path = dorado_bam_path
        self.mapped_bam = None
        self.dorado_bam = None

    def __enter__(self):
        self.mapped_bam = pysam.AlignmentFile(self.mapped_bam_path, "rb")
        self.dorado_bam = pysam.AlignmentFile(self.dorado_bam_path, "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.mapped_bam.close()
        self.dorado_bam.close()

    def _extract_methylation_data(self, read_id: str) -> List[Tuple[str, list]]:
        """Extracts methylation data for a given read ID from the Dorado BAM."""
        methylation_data = []
        for read in self.dorado_bam.fetch(until_eof=True):
            if read.query_name == read_id:
                # MM tag for modified bases, separated by bases
                meth_tag = read.get_tag('MM') if read.has_tag('MM') else ''
                # ML tag for methylation calls (probabilities)
                ml_tag = read.get_tag('ML') if read.has_tag('ML') else []
                methylation_data.append((meth_tag, ml_tag))
        return methylation_data

    def analyze(self) -> Dict[str, ReadData]:
        """Analyzes multiply-aligned reads, extracting methylation data and organizing by read ID."""
        reads_data = {}

        for read in self.mapped_bam:
            if not read.is_secondary:  # Consider only secondary alignments (multiply-mapped reads)
                continue

            read_id = read.query_name

            if read_id not in reads_data:
                reads_data[read_id] = ReadData(read_id)

            alignment_info = AlignmentInfo(
                reference_name=self.mapped_bam.get_reference_name(read.reference_id),
                reference_start=read.reference_start,
                reference_end=read.reference_end,
                query_alignment_start=read.query_alignment_start,
                query_alignment_end=read.query_alignment_end,
                cigar=read.cigarstring,
                mapping_quality=read.mapping_quality,
                flag=read.flag,
                methylation_data=self._extract_methylation_data(read_id)
            )
            reads_data[read_id].add_alignment(alignment_info)

        return reads_data

def main():
    mapped_bam_path = "path/to/your/mapped.bam"
    dorado_bam_path = "path/to/your/dorado.bam"

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