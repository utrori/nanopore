import pysam
import matplotlib.pyplot as plt
import pathlib as Path
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional

@dataclass
class AlignmentInfo:
    """Information about a specific alignment segment within a read."""
    reference_name: str
    reference_start: int
    reference_end: int
    query_alignment_start: int
    query_alignment_end: int
    length: int
    aligned_segment: pysam.AlignedSegment = None

@dataclass
class ReadData:
    """Information about a read including its alignments and read-specific attributes."""
    read_id: str
    cigar: str = None
    mapping_quality: int = None
    flag: int = None
    methylation_data: List[Tuple[int, int]] = field(default_factory=list)
    alignments: List[AlignmentInfo] = field(default_factory=list)
    has_strand_switch: bool = False  # Flag indicating whether strand switching is detected
    strand_switch_position: Optional[int] = None  # Changed from list to Optional[int]
    sequence: str = None  # Added read sequence
    read_length: int = None  # Added read length

    def add_alignment(self, alignment: AlignmentInfo):
        """Add an alignment to the read's alignments list."""
        self.alignments.append(alignment)
        
    def get_quality_scores(self):
        """Get quality scores from the primary alignment if available."""
        for alignment in self.alignments:
            if alignment.aligned_segment and not alignment.aligned_segment.is_secondary:
                return alignment.aligned_segment.query_qualities
        return None
    
    def analyze_strand_switching(self):
        """
        Analyze alignments to determine if the read exhibits strand switching.
        This method should be called after all alignments have been added.
        
        Returns:
            bool: True if strand switching is detected, False otherwise
        """
        # Reset strand switch data before analysis
        self.has_strand_switch = False
        self.strand_switch_position = None
        
        # We need at least 2 alignments to detect strand switching
        if len(self.alignments) < 2:
            return False
        
        # Sort alignments by query position for sequential analysis
        sorted_alignments = sorted(self.alignments, key=lambda a: a.query_alignment_start)
        
        # Analyze neighboring alignment pairs for potential strand switching
        for i in range(len(sorted_alignments) - 1):
            curr_alignment = sorted_alignments[i]
            next_alignment = sorted_alignments[i + 1]
            
            # Check if alignments are on different strands
            is_different_strand = curr_alignment.aligned_segment.is_reverse != next_alignment.aligned_segment.is_reverse
            
            # Calculate combined coverage of these two alignments
            curr_length = curr_alignment.query_alignment_end - curr_alignment.query_alignment_start
            next_length = next_alignment.query_alignment_end - next_alignment.query_alignment_start
            combined_coverage = (curr_length + next_length) / self.read_length if self.read_length else 0
            
            # If they're on different strands and cover more than 50% of the read
            if is_different_strand and combined_coverage > 0.5:
                self.has_strand_switch = True
                self.strand_switch_position = next_alignment.query_alignment_start
                break
        
        return self.has_strand_switch

class BamAnalyzer:
    # Mapped bam needs to be cleaned with samtools view -F 4
    def __init__(self, mapped_bam_path: str, dorado_bam_path: str):
        """
        Initialize analyzer with paths to mapped BAM and Dorado BAM files.
        
        Args:
            mapped_bam_path: Path to the BAM file with mapping information
            dorado_bam_path: Path to the BAM file with methylation data from Dorado
        """
        self.mapped_bam_path = mapped_bam_path
        self.dorado_bam_path = dorado_bam_path
        self.mapped_bam = None
        self.dorado_bam = None

    def __enter__(self):
        """Open BAM files when entering context."""
        self.dorado_bam = pysam.AlignmentFile(self.dorado_bam_path, "rb", check_sq=False)
        self.mapped_bam = pysam.AlignmentFile(self.mapped_bam_path, "rb", check_sq=False)
        self.methylation_data_by_read_id = self._load_methylation_data()  # Remove quality_data_by_read_id
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close BAM files when exiting context."""
        if self.mapped_bam:
            self.mapped_bam.close()
        if self.dorado_bam:
            self.dorado_bam.close()

    def _load_methylation_data(self) -> Dict[str, List[Tuple[int, int]]]:
        """
        Load methylation data from Dorado BAM.
        
        Returns:
            Dictionary mapping read IDs to methylation data
        """
        methylation_data_by_read_id = {}
        
        for read in self.dorado_bam:
            read_id = read.query_name
            if read.has_tag('MM') or read.has_tag('Ml'):
                modified_bases = read.modified_bases
                if ('C', 0, 'm') in modified_bases:
                    methylation_data_by_read_id[read_id] = read.modified_bases[('C', 0, 'm')]
        
        return methylation_data_by_read_id

    def analyze(self) -> Dict[str, ReadData]:
        """
        Analyze reads to extract mapping and methylation information.
        
        Returns:
            Dictionary mapping read IDs to ReadData objects
        """
        reads_data = {}

        for aligned_segment in self.mapped_bam:
            read_id = aligned_segment.query_name

            if read_id not in reads_data:
                # Create ReadData object with read-level attributes
                reads_data[read_id] = ReadData(
                    read_id=read_id,
                    cigar=aligned_segment.cigarstring,
                    mapping_quality=aligned_segment.mapping_quality,
                    flag=aligned_segment.flag,
                    methylation_data=self.methylation_data_by_read_id.get(read_id, []),
                    sequence=aligned_segment.query_sequence,  # Added sequence
                    read_length=aligned_segment.query_length  # Added read length
                )
            
            query_len = aligned_segment.query_length
            if aligned_segment.is_reverse:
                query_alignment_start = query_len - aligned_segment.query_alignment_end
                query_alignment_end = query_len - aligned_segment.query_alignment_start
            else:
                query_alignment_start = aligned_segment.query_alignment_start
                query_alignment_end = aligned_segment.query_alignment_end

            # Create alignment info
            alignment_info = AlignmentInfo(
                reference_name=self.mapped_bam.get_reference_name(aligned_segment.reference_id),
                reference_start=aligned_segment.reference_start,
                reference_end=aligned_segment.reference_end,
                query_alignment_start=query_alignment_start,
                query_alignment_end=query_alignment_end,
                length=aligned_segment.reference_end - aligned_segment.reference_start,
                aligned_segment=aligned_segment
            )
            reads_data[read_id].add_alignment(alignment_info)

        # After processing all alignments, analyze strand switching for each read
        for read_data in reads_data.values():
            if len(read_data.alignments) > 1:  # Only analyze reads with multiple alignments
                read_data.analyze_strand_switching()

        return reads_data

def methylation_per_coding():
    """
    Analyze methylation patterns in coding regions and create histogram.
    """
    dorado_bam_path = "test_files/R10.4.1_PSCA0047.bam"
    mapped_bam_path = "test_files/R10.4.1_PSCA0047/R10.4.1_PSCA0047.bam"
    with BamAnalyzer(mapped_bam_path, dorado_bam_path) as analyzer:
        reads_data = analyzer.analyze()
    
    ave_met_scores = []
    for read_id, read_data in reads_data.items():
        for alignment in read_data.alignments:
            met_list = read_data.methylation_data
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
    """Example usage of BamAnalyzer"""
    methylation_per_coding()
    quit()
    
    dorado_bam_path = "test_files/201020_47/calls_2025-02-05_T09-56-59.bam"
    mapped_bam_path = "test_files/201020_47_coding_mapped.bam"

    with BamAnalyzer(mapped_bam_path, dorado_bam_path) as analyzer:
        reads_data = analyzer.analyze()

    # Example of how to access and process the results
    for read_id, read_data in reads_data.items():
        print(f"Read ID: {read_id}")
        print(f"  CIGAR: {read_data.cigar}")
        print(f"  MAPQ: {read_data.mapping_quality}")
        print(f"  Flag: {read_data.flag}")
        print(f"  Sequence length: {read_data.read_length}")
        
        # Print the first 50 bases of the sequence (if available)
        if read_data.sequence:
            print(f"  Sequence (first 50bp): {read_data.sequence[:50]}...")
        
        # Get quality scores from the primary alignment if needed
        quality_scores = read_data.get_quality_scores()
        if quality_scores is not None:
            print(f"  Average quality: {np.mean(quality_scores):.2f}")
        
        if read_data.methylation_data:
            print(f"  Methylation data count: {len(read_data.methylation_data)}")
        else:
            print("  No methylation data found.")
            
        for i, alignment in enumerate(read_data.alignments):
            print(f"  Alignment {i+1}:")
            print(f"    Reference: {alignment.reference_name}:{alignment.reference_start}-{alignment.reference_end}")
            print(f"    Query Alignment: {alignment.query_alignment_start}-{alignment.query_alignment_end}")
            print(f"    Length: {alignment.length}")
        print("-" * 20)

    # Example of accessing strand switching information
    with BamAnalyzer(mapped_bam_path, dorado_bam_path) as analyzer:
        reads_data = analyzer.analyze()

    strand_switch_count = 0
    for read_id, read_data in reads_data.items():
        if read_data.has_strand_switch:
            strand_switch_count += 1
            print(f"Read {read_id} has strand switching at position: {read_data.strand_switch_position}")
    
    print(f"Total reads with strand switching: {strand_switch_count}/{len(reads_data)}")

if __name__ == "__main__":
    main()