import h5py
from pathlib import Path
import numpy as np
from typing import Dict, Tuple, Optional

class SingleFast5Reader:
    def __init__(self, filepath: str, bc_n: int = 0, methylation: str = 'cpg'):
        """Initialize Fast5 reader with file path."""
        self.filepath = Path(filepath)
        self.bc_n = bc_n
        self.methylation = methylation
        
    def extract_read_data(self) -> Dict:
        """Read single-read Fast5 file and extract key information."""
        try:
            with h5py.File(self.filepath, 'r') as f:
                # Navigate to the read group (assumes single read)
                """
                read_group = list(f['read'].keys())[0]
                read_path = f'/read/{read_group}'
                """
                
                # Extract sequence data
                basecall_group = f[f'/Analyses/Basecall_1D_00{self.bc_n}']
                fastq = basecall_group['BaseCalled_template/Fastq'][()].decode('utf-8')
                
                # Parse FASTQ format
                lines = fastq.split('\n')
                sequence = lines[1]
                quality_string = lines[3]
                phred_scores = [ord(c) - 33 for c in quality_string]
                
                # Get methylation calls if available
                if self.methylation:
                    methylation_data = self._get_methylation_calls(f, self.methylation)
                
                return {
                    'read_id': self.filepath.stem,
                    'sequence': sequence,
                    'quality_string': quality_string,
                    'quality_scores': phred_scores,
                    'methylation': methylation_data,
                }
                
        except Exception as e:
            raise ValueError(f"Error reading Fast5 file: {str(e)}")
    
    def _get_methylation_calls(self, f: h5py.File, methylation: str) -> Optional[list]:
        """Extract methylation calls if available."""
        try:
            basecall_group = f[f'/Analyses/Basecall_1D_00{self.bc_n}']
            mod_base_table = basecall_group['BaseCalled_template/ModBaseProbs'][:]
            if methylation == 'cpg':
                return mod_base_table[:,2]
            if methylation == '6ma':
                return mod_base_table[:,1]
        except:
            pass
        return None

def main():
    # Example usage
    reader = SingleFast5Reader("/mnt/data/nanopore_cas9_bc/221124_Y8585fN_CR_bc/workspace/0/00128cf9-6898-4859-a5dd-b94bb8b4d324.fast5")
    data = reader.extract_read_data()
    print(f"Read ID: {data['read_id']}")
    print(f"Sequence length: {len(data['sequence'])}")
    print(f"Quality scores length: {len(data['quality_scores'])}")
    if data['methylation'] is not None:
        print(len(data['methylation']))

if __name__ == "__main__":
    main()