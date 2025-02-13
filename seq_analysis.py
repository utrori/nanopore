from fast5_loader import SingleFast5Reader
from dataclasses import dataclass
import utilities
import numpy as np
import config

@dataclass
class SplitAnalyzedRead:
    ref: str
    read_id: str
    methylation_data: list
    ave_quality: float
    length: int
    split_mapping_res: np.array

@dataclass
class MinimapAnalyzedRead:
    read_id: str
    methylation_data: list
    ave_quality: float
    length: int
    mapped_regions: list
    cigar_converted: list

class ReadAnalyzer:
    """Analyzes a single read from a FAST5 file."""
    def __init__(self, reader: object, ref: str):
        """Initialize read analyzer with FAST5 file path."""
        self.ref = ref
        self.loadeddata = reader.extract_read_data()
        self.analyzed_read = None

    def split_alignment(self) -> SplitAnalyzedRead:
        """Split the read into multiple parts and align each part."""
        seq = self.loadeddata['sequence']
        qual = self.loadeddata['quality_string']
        split_mapping_res = utilities.split_mapping_and_sam_analysis(config.SPLIT_LENGTH, self.loadeddata['read_id'], seq, qual, self.ref)    
        methylation_data = self.loadeddata['methylation']
        ave_quality = np.mean(self.loadeddata['quality_scores'])
        length = len(self.loadeddata['sequence'])
        self.analyzed_read = SplitAnalyzedRead(
            ref = self.ref,
            read_id=self.loadeddata['read_id'],
            methylation_data=methylation_data,
            ave_quality=ave_quality,
            length=length,
            split_mapping_res=split_mapping_res
        )
        return self.analyzed_read
    
    def minimap2_alignment(self) -> MinimapAnalyzedRead:
        """Align the read using minimap2."""
        seq = self.loadeddata['sequence']
        qual = self.loadeddata['quality_string']
        mapped_regions, cigar_converted = utilities.minimap_mapping_and_sam_analysis(self.loadeddata['read_id'], seq, qual, self.ref)
        methylation_data = self.loadeddata['methylation']
        ave_quality = np.mean(self.loadeddata['quality_scores'])
        length = len(self.loadeddata['sequence'])
        self.analyzed_read = MinimapAnalyzedRead(
            read_id=self.loadeddata['read_id'],
            methylation_data=methylation_data,
            ave_quality=ave_quality,
            length=length,
            mapped_regions=mapped_regions,
            cigar_converted=cigar_converted
        )
        return self.analyzed_read