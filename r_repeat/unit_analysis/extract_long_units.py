"""
Utility script to extract sequences for unusually long r-repeat units.
"""

import argparse
import json
import pandas as pd
from pathlib import Path
import pysam
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def extract_unit_sequences(long_units_csv, bam_path, output_dir=None):
    """
    Extract sequence data for identified long r-repeat units.
    
    Parameters:
        long_units_csv: CSV file containing long unit information
        bam_path: Path to BAM file containing original reads
        output_dir: Directory for output files
        
    Returns:
        str: Path to output FASTA file
    """
    # Read the CSV file
    df = pd.read_csv(long_units_csv)
    
    if df.empty:
        logger.error("No long units found in the CSV file")
        return None
    
    if output_dir is None:
        output_dir = Path("results/long_units")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_fasta = output_dir / "long_unit_sequences.fa"
    
    extracted_count = 0
    
    with open(output_fasta, 'w') as f_out:
        # Open the BAM file
        with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
            # Process each long unit
            for _, row in df.iterrows():
                read_id = row['read_id']
                start = int(row['start'])
                end = int(row['end'])
                unit_length = int(row['unit_length'])
                region_idx = int(row['region_index'])
                unit_idx = int(row['unit_index'])
                
                # Find the read
                found = False
                for read in bam.fetch(until_eof=True):
                    if read.query_name == read_id:
                        if read.query_sequence:
                            # Extract the unit sequence
                            sequence = read.query_sequence[start:end]
                            
                            # Write to FASTA
                            f_out.write(f">{read_id}_r{region_idx+1}_u{unit_idx+1}_{unit_length}bp\n")
                            
                            # Write sequence with line breaks every 80 characters
                            for i in range(0, len(sequence), 80):
                                f_out.write(f"{sequence[i:i+80]}\n")
                            
                            extracted_count += 1
                            found = True
                        break
                
                if not found:
                    logger.warning(f"Read {read_id} not found in BAM file")
    
    logger.info(f"Extracted {extracted_count}/{len(df)} long unit sequences to {output_fasta}")
    return str(output_fasta)

def test_extract_long_units(bam_file=None, r_repeat_regions_json=None, threshold=1000, output_dir=None):
    """
    Test function to find and extract long r-repeat units without requiring argument parsing.
    
    This function directly uses the unit length information that's already available
    in the r-repeat regions data structure to identify and extract long units.
    
    Parameters:
        bam_file: Path to the BAM file containing read sequences
        r_repeat_regions_json: Path to JSON file with r-repeat region data
        threshold: Minimum length to consider a unit "long" (default: 1000bp)
        output_dir: Directory for output files (default: "results/long_units")
        
    Returns:
        tuple: (CSV file with long unit info, FASTA file with long unit sequences)
    """
    import sys
    
    # Add parent directory to path to import the unit analysis module if needed
    parent_dir = str(Path(__file__).parent)
    if parent_dir not in sys.path:
        sys.path.append(parent_dir)
    
    logger.info("Starting long r-repeat unit extraction test")
    
    # Set default paths if not provided
    if bam_file is None:
        # Try to find a BAM file in common locations
        potential_paths = [
            "/home/owner/nanopore_rDNA/data/new_r_repeat_analysis/unit_analysis/HG00438_1_motif_aligned.bam",
            "/media/owner/bbeedd59-668c-4d67-a65b-442d3272c3bd/hpgp/dorado_bc/HG00733.bam"
        ]
        
        for path in potential_paths:
            if Path(path).exists():
                bam_file = path
                break
                
        if bam_file is None:
            logger.error("No BAM file found. Please specify a BAM file path.")
            return None, None
    
    if r_repeat_regions_json is None:
        # Try to find a regions JSON file in common locations
        bam_name = Path(bam_file).stem
        potential_paths = [
            "/home/owner/nanopore_rDNA/data/new_r_repeat_analysis/r_repeat_regions.json",
            f"r_repeat_results/{bam_name}/r_repeat_regions.json",
            f"results/{bam_name}/r_repeat_regions.json",
            "r_repeat_results/latest/r_repeat_regions.json"
        ]
        
        for path in potential_paths:
            if Path(path).exists():
                r_repeat_regions_json = path
                break
                
        if r_repeat_regions_json is None:
            logger.error("No r-repeat regions JSON file found. Please specify a JSON file path.")
            return None, None
    
    if output_dir is None:
        output_dir = Path("data/new_r_repeat_analysis/long_units")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load r-repeat regions
    logger.info(f"Loading r-repeat regions from {r_repeat_regions_json}")
    with open(r_repeat_regions_json, 'r') as f:
        r_repeat_regions = json.load(f)
    
    # Find long units directly from the unit length information or using the has_long_unit flag
    long_units = []
    
    # Process the r-repeat regions to find long units
    for read_id, regions in r_repeat_regions.items():
        # Handle both single region and multiple regions per read
        if isinstance(regions, list):
            regions_list = regions
        else:
            regions_list = [regions]
            
        for region_idx, r_repeat_region in enumerate(regions_list):
            # First check if the region is already marked as having long units
            if r_repeat_region.get('has_long_unit', False):
                # If unit_lengths is available, find which ones are long
                if 'unit_lengths' in r_repeat_region:
                    unit_lengths = r_repeat_region['unit_lengths']
                    
                    # Find long units
                    for unit_idx, unit_length in enumerate(unit_lengths):
                        if unit_length > threshold:
                            # For comprehensive JSON format, we already have unit lengths 
                            # and can estimate unit positions more accurately
                            
                            # Calculate start position based on unit index
                            if unit_idx == 0:
                                # First unit starts at beginning of region
                                unit_start = r_repeat_region['start']
                            else:
                                # Use start + sum of lengths of previous units
                                unit_start = r_repeat_region['start'] + sum(unit_lengths[:unit_idx])
                            
                            # End position is start + this unit's length
                            unit_end = unit_start + unit_length
                            
                            # Add to long units list
                            long_units.append({
                                'read_id': read_id,
                                'region_index': region_idx,
                                'unit_index': unit_idx,
                                'unit_length': unit_length,
                                'start': unit_start,
                                'end': unit_end
                            })
    
    # Sort units by length (descending)
    long_units.sort(key=lambda x: x['unit_length'], reverse=True)
    
    logger.info(f"Found {len(long_units)} units longer than {threshold}bp")
    
    # Print the top 10 longest units
    if long_units:
        logger.info("Top 10 longest units:")
        for i, unit in enumerate(long_units[:10]):
            logger.info(f"{i+1}. Read {unit['read_id']} (region {unit['region_index']}, "
                      f"unit {unit['unit_index']}): {unit['unit_length']}bp")
    
    # Output to CSV
    csv_file = output_dir / f"long_units_over_{threshold}bp.csv"
    with open(csv_file, 'w') as f:
        # Write header
        f.write("read_id,region_index,unit_index,unit_length,start,end\n")
        
        # Write data
        for unit in long_units:
            f.write(f"{unit['read_id']},{unit['region_index']},{unit['unit_index']},"
                   f"{unit['unit_length']},{unit['start']},{unit['end']}\n")
    
    logger.info(f"Test completed. Found {len(long_units)} long units.")
    logger.info(f"Results saved to {csv_file}")
    
    return csv_file

def main():
    test_extract_long_units()
    quit()
    parser = argparse.ArgumentParser(description="Extract sequences for long r-repeat units")
    parser.add_argument("--units-csv", required=True, help="CSV file containing long unit information")
    parser.add_argument("--bam", required=True, help="BAM file containing original reads")
    parser.add_argument("--output-dir", help="Directory for output files")
    
    args = parser.parse_args()
    
    extract_unit_sequences(args.units_csv, args.bam, args.output_dir)

if __name__ == "__main__":
    # Test function that can be called directly without argparse
    import sys
    if len(sys.argv) == 1:
        # No arguments provided, run the test function
        test_extract_long_units()
    else:
        # Arguments provided, use the regular argparse flow
        main()
