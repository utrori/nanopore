"""
R-repeat region detection module.

This module detects r-repeat regions in Oxford Nanopore sequencing reads by:
1. Creating FASTQ files from BAM inputs
2. Using minimap2 to align reads against r-repeat references
3. Identifying reads with r-repeat regions and their boundaries
"""

from pathlib import Path
import os
import json
import pysam
import subprocess
import logging
import tempfile
import sys

# Ensure parent directory is in path for imports
parent_dir = str(Path(__file__).parent.parent.parent)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import numpy as np  # Add numpy import
import shutil  # Add shutil import
import argparse  # Add argparse import

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler("logs/r_repeat_detection.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def create_flanking_references(output_dir=None):
    """
    Create reference files for the 5' and 3' flanking regions of r-repeats.
    
    Parameters:
        output_dir: Directory for output files
        
    Returns:
        tuple: Paths to 5' and 3' flank reference files
    """
    if output_dir is None:
        output_dir = Path("references/flanking_regions")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get reference sequence
    human_rDNA_ref = Path("references/human_rDNA.fa")
    seq = ""
    with open(human_rDNA_ref, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq += line.strip()
    
    # Extract the flanking sequences (based on known coordinates)
    flank_5_seq = seq[12993:13493]  # 500bp upstream of r-repeat
    flank_3_seq = seq[15584:16084]  # 500bp downstream of r-repeat
    
    # Create reference files
    flank_5_file = output_dir / "r_repeat_5_flank.fa"
    with open(flank_5_file, 'w') as f:
        f.write(f">r_repeat_5_flank\n{flank_5_seq}\n")
    
    flank_3_file = output_dir / "r_repeat_3_flank.fa"
    with open(flank_3_file, 'w') as f:
        f.write(f">r_repeat_3_flank\n{flank_3_seq}\n")
    
    logger.info(f"Created flanking region references in {output_dir}")
    return str(flank_5_file), str(flank_3_file)

def minimap2_mapping(query_file, reference_file, out_bam=None, threads=4):
    """
    Perform alignment using minimap2 and convert to BAM.
    
    Parameters:
        query_file: Path to query sequences (FASTQ or FASTA)
        reference_file: Path to reference sequence
        out_bam: Output BAM file path (default: based on input names)
        threads: Number of threads to use
        
    Returns:
        str: Path to output BAM file
    """
    if out_bam is None:
        query_name = Path(query_file).stem
        ref_name = Path(reference_file).stem
        out_bam = f"temp_files/{query_name}_to_{ref_name}.bam"
    
    # Ensure output directory exists
    Path(out_bam).parent.mkdir(parents=True, exist_ok=True)
    
    # Run minimap2 and pipe to samtools
    cmd = [
        "minimap2", "-ax", "map-ont", 
        "-Y", "-t", str(threads),
        str(reference_file), str(query_file)  # Convert Path objects to strings
    ]
    
    with open(out_bam, 'wb') as bam_file:
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ["samtools", "view", "-b", "-"],
            stdin=p1.stdout, stdout=subprocess.PIPE
        )
        p3 = subprocess.Popen(
            ["samtools", "sort", "-o", str(out_bam)],  # Convert Path to string
            stdin=p2.stdout
        )
        p1.stdout.close()
        p2.stdout.close()
        p3.wait()
    
    # Index the BAM file - convert Path to string
    pysam.index(str(out_bam))
    
    logger.info(f"Aligned {query_file} to {reference_file}, output: {out_bam}")
    return str(out_bam)  # Return string path instead of Path object

def align_reads_to_flanking_regions(fastq_path, output_dir=None, reuse_existing=False):
    """
    Align reads to r-repeat flanking regions to identify r-repeat boundaries.
    
    Parameters:
        fastq_path: Path to FASTQ file
        output_dir: Directory for output files
        reuse_existing: Whether to reuse existing alignments if available
        
    Returns:
        tuple: Paths to 5' and 3' flank alignment BAM files
    """
    if output_dir is None:
        output_dir = Path("temp_files/flanking_alignments")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get references
    flank_5_ref, flank_3_ref = create_flanking_references()
    
    # Define output BAM paths
    flank_5_bam = output_dir / f"{Path(fastq_path).stem}_5flank.bam"
    flank_3_bam = output_dir / f"{Path(fastq_path).stem}_3flank.bam"
    
    # Check if we can reuse existing alignments
    if reuse_existing and flank_5_bam.exists() and flank_3_bam.exists() and Path(f"{flank_5_bam}.bai").exists() and Path(f"{flank_3_bam}.bai").exists():
        logger.info("Using existing flanking region alignments")
        return str(flank_5_bam), str(flank_3_bam)
    
    # Align to 5' flank
    minimap2_mapping(fastq_path, flank_5_ref, out_bam=flank_5_bam)
    
    # Align to 3' flank
    minimap2_mapping(fastq_path, flank_3_ref, out_bam=flank_3_bam)
    
    logger.info(f"Aligned reads to flanking regions: {flank_5_bam} and {flank_3_bam}")
    return str(flank_5_bam), str(flank_3_bam)

def identify_r_repeat_regions(flank_5_bam, flank_3_bam, allow_multiple=True):
    """
    Identify r-repeat regions in nanopore reads based on flanking region alignments.
    
    This function analyzes alignments to the 5' and 3' flanking regions to determine
    the boundaries of r-repeats within reads, accounting for strand orientation and
    filtering out implausible r-repeat lengths.
    
    Parameters:
        flank_5_bam: Path to BAM file with 5' flank alignments
        flank_3_bam: Path to BAM file with 3' flank alignments
        allow_multiple: Whether to allow multiple r-repeat regions per read
        
    Returns:
        dict: Dictionary mapping read IDs to r-repeat regions (either a single dict or list of dicts)
    """
    # Process 5' flank alignments - store all alignments
    read_id_to_5flank = {}
    with pysam.AlignmentFile(flank_5_bam, "rb") as bam:
        for aln in bam:
            if aln.is_unmapped or aln.mapping_quality < 20 or np.mean(aln.query_qualities) < 16:
                continue
            
            # Record this 5' flank alignment
            if allow_multiple:
                if aln.query_name not in read_id_to_5flank:
                    read_id_to_5flank[aln.query_name] = []
                
                read_id_to_5flank[aln.query_name].append({
                    'r_start': aln.query_alignment_end,
                    'is_reverse': aln.is_reverse,
                    'query_length': aln.query_length,
                    'mapping_quality': aln.mapping_quality,
                    'alignment_score': aln.get_tag('AS') if aln.has_tag('AS') else 0
                })
            else:
                # Original behavior - keep only the best alignment
                if (aln.query_name not in read_id_to_5flank or 
                    aln.mapping_quality > read_id_to_5flank[aln.query_name]['mapping_quality']):
                    read_id_to_5flank[aln.query_name] = {
                        'r_start': aln.query_alignment_end,
                        'is_reverse': aln.is_reverse,
                        'query_length': aln.query_length,
                        'mapping_quality': aln.mapping_quality
                    }
    
    # Process 3' flank alignments - store all alignments
    read_id_to_3flank = {}
    with pysam.AlignmentFile(flank_3_bam, "rb") as bam:
        for aln in bam:
            if aln.is_unmapped or aln.mapping_quality < 20:
                continue
            
            # Record this 3' flank alignment
            if allow_multiple:
                if aln.query_name not in read_id_to_3flank:
                    read_id_to_3flank[aln.query_name] = []
                
                read_id_to_3flank[aln.query_name].append({
                    'r_end': aln.query_alignment_start,
                    'is_reverse': aln.is_reverse,
                    'query_length': aln.query_length,
                    'mapping_quality': aln.mapping_quality,
                    'alignment_score': aln.get_tag('AS') if aln.has_tag('AS') else 0
                })
            else:
                # Original behavior - keep only the best alignment
                if (aln.query_name not in read_id_to_3flank or 
                    aln.mapping_quality > read_id_to_3flank[aln.query_name]['mapping_quality']):
                    read_id_to_3flank[aln.query_name] = {
                        'r_end': aln.query_alignment_start,
                        'is_reverse': aln.is_reverse,
                        'query_length': aln.query_length,
                        'mapping_quality': aln.mapping_quality
                    }
    
    # Find r-repeat regions
    read_id_to_r_repeat = {}
    
    # Process different based on whether we allow multiple regions per read
    if allow_multiple:
        # Look for multiple r-repeats per read
        for read_id in set(read_id_to_5flank.keys()) & set(read_id_to_3flank.keys()):
            # Initialize list for this read
            read_id_to_r_repeat[read_id] = []
            
            # Try all combinations of 5' and 3' flank alignments
            for flank_5 in read_id_to_5flank[read_id]:
                for flank_3 in read_id_to_3flank[read_id]:
                    # Both flanks must be on the same strand
                    if flank_5['is_reverse'] != flank_3['is_reverse']:
                        continue
                    
                    # Both flanks should have comparable query_length
                    if abs(flank_5['query_length'] - flank_3['query_length']) > 100:
                        continue
                    
                    is_reverse = flank_5['is_reverse']
                    read_length = flank_5['query_length']
                    
                    # Get raw boundaries
                    r_start = flank_5['r_start']
                    r_end = flank_3['r_end']
                    
                    """
                    # Adjust positions for reverse strand
                    if is_reverse:
                        old_start, old_end = r_start, r_end
                        r_start = read_length - old_end
                        r_end = read_length - old_start
                    """
                    
                    # Check valid order and minimum separation distance
                    if r_start >= r_end or r_end - r_start < 500:
                        continue
                    
                    # Calculate r-repeat length
                    r_length = r_end - r_start
                    
                    # Only include reasonable r-repeat lengths
                    if 500 <= r_length <= 5000:
                        # Create r-repeat region entry
                        region = {
                            'start': r_start,
                            'end': r_end,
                            'length': r_length,
                            'is_reverse': is_reverse,
                            'flank_5_quality': flank_5['mapping_quality'],
                            'flank_3_quality': flank_3['mapping_quality'],
                            'combined_quality': flank_5['mapping_quality'] + flank_3['mapping_quality']
                        }
                        
                        # Add alignment scores if available
                        if 'alignment_score' in flank_5 and 'alignment_score' in flank_3:
                            region['combined_score'] = flank_5['alignment_score'] + flank_3['alignment_score']
                        
                        # Add to read's regions
                        read_id_to_r_repeat[read_id].append(region)
            
            # Sort regions by start position
            if read_id_to_r_repeat[read_id]:
                # Sort by start_position
                read_id_to_r_repeat[read_id].sort(
                    key=lambda x: x['start']) 
            # Add index to each region
            for i, region in enumerate(read_id_to_r_repeat[read_id]):
                region['index'] = i
                
            """
                # Remove overlapping regions (keep higher quality ones)
                non_overlapping = []
                for region in read_id_to_r_repeat[read_id]:
                    # Check if this region overlaps with any already selected regions
                    overlaps = False
                    for selected in non_overlapping:
                        # Check for overlap
                        if (region['start'] < selected['end'] and region['end'] > selected['start']):
                            overlaps = True
                            break
                    
                    # Add if no overlap
                    if not overlaps:
                        non_overlapping.append(region)
                
                # Update with non-overlapping set
                read_id_to_r_repeat[read_id] = non_overlapping
            """
            
            # Remove reads with no valid regions
            if not read_id_to_r_repeat[read_id]:
                del read_id_to_r_repeat[read_id]
    else:
        # Original behavior - find a single r-repeat per read
        for read_id in set(read_id_to_5flank.keys()) & set(read_id_to_3flank.keys()):
            # Both flanks should be on the same strand
            if read_id_to_5flank[read_id]['is_reverse'] != read_id_to_3flank[read_id]['is_reverse']:
                continue
            
            # Get r-repeat boundaries
            r_start = read_id_to_5flank[read_id]['r_start']
            r_end = read_id_to_3flank[read_id]['r_end']
            is_reverse = read_id_to_5flank[read_id]['is_reverse']
            
            # Adjust positions for reverse strand directly here
            if is_reverse:
                read_length = read_id_to_5flank[read_id]['query_length']
                # Flip the positions for reverse strand reads
                old_start, old_end = r_start, r_end
                r_start = read_length - old_end
                r_end = read_length - old_start
            
            # Check that boundaries make sense (r_start should be before r_end)
            if r_start >= r_end:
                continue
                
            # Calculate r-repeat length
            r_length = r_end - r_start
            
            # Only include reasonable r-repeat lengths
            if 500 <= r_length <= 5000:
                read_id_to_r_repeat[read_id] = {
                    'start': r_start,
                    'end': r_end,
                    'length': r_length,
                    'is_reverse': is_reverse
                }
    
    # Log statistics
    total_regions = 0
    lengths = []
    
    if allow_multiple:
        # Count total regions across all reads
        for regions in read_id_to_r_repeat.values():
            total_regions += len(regions)
            lengths.extend([region['length'] for region in regions])
    else:
        total_regions = len(read_id_to_r_repeat)
        lengths = [data['length'] for data in read_id_to_r_repeat.values()]
    
    if lengths:
        logger.info(f"Identified {total_regions} r-repeat regions in {len(read_id_to_r_repeat)} reads")
        logger.info(f"Length statistics: min={min(lengths)}, median={np.median(lengths):.1f}, max={max(lengths)}")
        
        # Calculate reads with multiple regions
        if allow_multiple:
            multi_region_reads = sum(1 for regions in read_id_to_r_repeat.values() if len(regions) > 1)
            logger.info(f"Found {multi_region_reads} reads ({multi_region_reads/len(read_id_to_r_repeat):.1%}) with multiple r-repeat regions")
    else:
        logger.warning("No r-repeat regions found")
    
    return read_id_to_r_repeat

def extract_r_repeat_sequences(read_id_to_r_repeat, fastq_path, output_file=None):
    """
    Extract r-repeat sequences from reads based on identified coordinates.
    
    Parameters:
        read_id_to_r_repeat: Dictionary mapping read IDs to r-repeat regions
        fastq_path: Path to FASTQ file containing reads
        output_file: Path for output FASTA file
        
    Returns:
        str: Path to the FASTA file with r-repeat sequences
    """
    if output_file is None:
        sample_name = Path(fastq_path).stem
        output_file = f"temp_files/r_repeat_sequences/{sample_name}_r_repeats.fa"
    
    # Ensure output directory exists
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    # Extract sequences - using different approaches depending on file type
    extracted_count = 0
    
    # Determine if we have multiple regions per read
    has_multiple_regions = any(isinstance(regions, list) for regions in read_id_to_r_repeat.values())
    
    with open(output_file, 'w') as f_out:
        # Check if input is FASTQ format
        if str(fastq_path).lower().endswith(('.fastq', '.fq')):
            # Parse FASTQ file manually
            with open(fastq_path, 'r') as f_in:
                line_count = 0
                current_read_id = None
                current_sequence = None
                
                for line in f_in:
                    line_count += 1
                    line_type = line_count % 4
                    
                    if line_type == 1:  # Header line with read ID
                        # Extract read ID without @ prefix and comments
                        current_read_id = line.strip().split()[0][1:]
                    elif line_type == 2:  # Sequence line
                        current_sequence = line.strip()
                    elif line_type == 0:  # End of read entry
                        # Process the read if it's in our list
                        if current_read_id in read_id_to_r_repeat:
                            # Handle either single or multiple regions
                            if has_multiple_regions:
                                regions = read_id_to_r_repeat[current_read_id]
                                for i, region in enumerate(regions):
                                    r_start = region['start']
                                    r_end = region['end']
                                    is_reverse = region['is_reverse']
                                    
                                    # Extract sequence
                                    if not is_reverse:
                                        sequence = current_sequence[r_start:r_end]
                                    else:
                                        sequence = reverse_complement(current_sequence)[r_start:r_end]
                                    
                                    
                                    # Write to FASTA with region index
                                    f_out.write(f">{current_read_id}_r{i+1}\n{sequence}\n")
                                    extracted_count += 1
                            else:
                                # Original behavior for single region
                                coords = read_id_to_r_repeat[current_read_id]
                                r_start = coords['start']
                                r_end = coords['end']
                                
                                # Extract sequence
                                sequence = current_sequence[r_start:r_end]
                                
                                # Write to FASTA
                                f_out.write(f">{current_read_id}\n{sequence}\n")
                                extracted_count += 1
        else:
            # Assume it's a BAM/CRAM file
            try:
                with pysam.AlignmentFile(fastq_path, "r", check_sq=False) as bam:
                    for read in bam:
                        read_id = read.query_name
                        if read_id not in read_id_to_r_repeat:
                            continue
                        
                        # Handle either single or multiple regions
                        if has_multiple_regions:
                            regions = read_id_to_r_repeat[read_id]
                            for i, region in enumerate(regions):
                                r_start = region['start']
                                r_end = region['end']
                                
                                # Extract sequence - check if read has sequence data
                                if not read.query_sequence:
                                    continue
                                
                                # Extract sequence using coordinates that are already adjusted for strand
                                sequence = read.query_sequence[r_start:r_end]
                                
                                # Write to FASTA with region index
                                f_out.write(f">{read_id}_r{i+1}\n{sequence}\n")
                                extracted_count += 1
                        else:
                            # Original behavior for single region
                            coords = read_id_to_r_repeat[read_id]
                            r_start = coords['start']
                            r_end = coords['end']
                            
                            # Extract sequence - check if read has sequence data
                            if not read.query_sequence:
                                continue
                            
                            # Extract sequence using coordinates that are already adjusted for strand
                            sequence = read.query_sequence[r_start:r_end]
                            
                            # Write to FASTA
                            f_out.write(f">{read_id}\n{sequence}\n")
                            extracted_count += 1
            except Exception as e:
                logger.error(f"Error processing file {fastq_path}: {str(e)}")
                logger.error("Please ensure the input file is a valid FASTQ or BAM file.")
                raise
    
    # Log results
    if has_multiple_regions:
        total_regions = sum(len(regions) for regions in read_id_to_r_repeat.values())
        logger.info(f"Extracted {extracted_count}/{total_regions} r-repeat sequences from {len(read_id_to_r_repeat)} reads")
    else:
        logger.info(f"Extracted {extracted_count}/{len(read_id_to_r_repeat)} r-repeat sequences")
    
    return output_file

def reverse_complement(sequence):
    """
    Get the reverse complement of a DNA sequence.
    
    Parameters:
        sequence: DNA sequence string
        
    Returns:
        str: Reverse complement sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'N': 'N', '.': '.', '-': '-', '*': '*'}
    return ''.join(complement.get(base.upper(), 'N') for base in reversed(sequence))

def detect_r_repeats(fastq_file, output_dir=None, reuse_existing=False, allow_multiple=True, temp_dir=None):
    """
    Main function to detect r-repeat regions in a FASTQ file.
    
    Parameters:
        fastq_file: Path to FASTQ file with reads
        output_dir: Directory for output files
        reuse_existing: Whether to reuse existing alignments if available
        allow_multiple: Whether to allow multiple r-repeat regions per read
        temp_dir: Directory for temporary files
        
    Returns:
        tuple: (dict of r-repeat regions, path to extracted sequences FASTA)
    """
    if output_dir is None:
        output_dir = Path("r_repeat_results") / Path(fastq_file).stem
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create data directory for serialized data
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)
    
    # Create temp directory if not provided
    if temp_dir is None:
        temp_dir = output_dir / "temp" / "detection"
    else:
        temp_dir = Path(temp_dir) / "detection"
        
    temp_dir.mkdir(exist_ok=True, parents=True)
    
    # Align reads to flanking regions
    flank_5_bam, flank_3_bam = align_reads_to_flanking_regions(
        fastq_file, 
        output_dir=temp_dir, 
        reuse_existing=reuse_existing
    )
    
    # Identify r-repeat regions with option for multiple regions per read
    r_repeat_regions = identify_r_repeat_regions(flank_5_bam, flank_3_bam, allow_multiple=allow_multiple)
    
    # Save r-repeat regions to file
    regions_file = data_dir / "r_repeat_regions.json"
    with open(regions_file, 'w') as f:
        # Convert r-repeat regions to a serializable format
        serializable_regions = {
            read_id: [
                {
                    'start': int(region['start']),
                    'end': int(region['end']),
                    'length': int(region['length']),
                    'is_reverse': bool(region['is_reverse'])
                }
                for region in regions
            ] if isinstance(regions, list) else {
                'start': int(regions['start']),
                'end': int(regions['end']),
                'length': int(regions['length']),
                'is_reverse': bool(regions['is_reverse'])
            }
            for read_id, regions in r_repeat_regions.items()
        }
        json.dump(serializable_regions, f, indent=2)
    
    # Extract r-repeat sequences
    sequences_file = data_dir / "r_repeat_sequences.fa"
    extract_r_repeat_sequences(r_repeat_regions, fastq_file, sequences_file)
    
    # Save a copy to the main output directory for easy access
    output_sequences_file = output_dir / "r_repeat_sequences.fa"
    shutil.copy(sequences_file, output_sequences_file)
    
    logger.info(f"R-repeat detection completed. Results saved to {output_dir}")
    return r_repeat_regions, str(output_sequences_file)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Detect r-repeat regions in nanopore reads")
    parser.add_argument("fastq", help="Path to FASTQ file containing reads")
    parser.add_argument("--output", "-o", help="Output directory for results")
    parser.add_argument("--reuse", "-r", action="store_true", help="Reuse existing alignments if available")
    parser.add_argument("--allow-multiple", "-m", action="store_true", help="Allow multiple r-repeat regions per read")
    
    args = parser.parse_args()
    
    detect_r_repeats(args.fastq, args.output, args.reuse, args.allow_multiple)
