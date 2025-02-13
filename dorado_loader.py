from pathlib import Path
import pysam
import numpy as np
from typing import Dict, Tuple, Optional

import pysam
from pathlib import Path
import numpy as np
from typing import Dict, Tuple, Optional, Iterator

class BamAlignmentWrapper:
    """Wraps a pysam.AlignedSegment object to mimic a SingleFast5Reader."""

    def __init__(self, alignment: pysam.AlignedSegment):
        """Initialize with a pysam.AlignedSegment object."""
        self.alignment = alignment

    def extract_read_data(self) -> Dict:
        """Mimics SingleFast5Reader.extract_read_data() using data from the AlignedSegment."""
        sequence = self.alignment.query_sequence
        quality_string = "".join([chr(q + 33) for q in self.alignment.query_qualities])
        phred_scores = self.alignment.query_qualities
        methylation_data = self.alignment.modified_bases[('C', 0, 'm')] if ('C', 0, 'm') in self.alignment.modified_bases else None

        return {
            'read_id': self.alignment.query_name,
            'sequence': sequence,
            'quality_string': quality_string,
            'quality_scores': phred_scores,
            'methylation': methylation_data,
        }


class BamToSingleReadReader:
    def __init__(self, bam_path: str):
        """
        Initializes a BAM reader that yields BamAlignmentWrapper objects.

        Args:
            bam_path: Path to the BAM file.
            bc_n: Basecall number (used for methylation data).
            methylation: Type of methylation to extract ('cpg' or '6ma').
        """
        self.bam_path = Path(bam_path)
        self.bamfile = None

    def __enter__(self):
        """Opens the BAM file for reading when used in a `with` statement."""
        self.bamfile = pysam.AlignmentFile(self.bam_path, "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Closes the BAM file."""
        if self.bamfile:
            self.bamfile.close()

    def __iter__(self) -> Iterator[BamAlignmentWrapper]:
        """Iterates through the reads in the BAM file, yielding a BamAlignmentWrapper object for each read."""
        if self.bamfile is None:
            raise ValueError("BAM file is not open. Use this class with a 'with' statement.")

        for read in self.bamfile:
            yield BamAlignmentWrapper(read)

def main():
    # Example Usage
    bam_path = "test_files/your_bam_file.bam"  # Replace with your BAM file path

    with BamToSingleReadReader(bam_path) as bam_reader:
        for single_read in bam_reader:
            data = single_read.extract_read_data()
            print(f"Read ID: {data['read_id']}")
            print(f"Sequence length: {len(data['sequence'])}")
            print(f"Quality scores length: {len(data['quality_scores'])}")

            if data['methylation']:
                print(f"Methylation calls: {len(data['methylation'])}")
            else:
                print("No methylation data found for this read.")

            print("-" * 20)

if __name__ == "__main__":
    main()