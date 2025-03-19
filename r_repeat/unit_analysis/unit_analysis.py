"""
R-repeat unit structure analysis module - New Implementation.

This module analyzes the internal structure of r-repeat regions by identifying
the repetitive units composed of:
1. A conserved ~400bp motif
2. A variable CT dinucleotide repeat region
"""

import logging
import subprocess
from pathlib import Path
import pysam
import numpy as np
import json  # Add json import
import matplotlib.pyplot as plt  # Add matplotlib import
import seaborn as sns  # Add seaborn import

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler("logs/r_repeat_unit_analysis_new.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

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

def create_motif_reference(output_dir=None):
    """
    Create reference file for the conserved ~400bp motif in r-repeats.
    
    Parameters:
        output_dir: Directory for output file
        
    Returns:
        str: Path to the motif reference file
    """
    if output_dir is None:
        output_dir = Path("references/r_repeat_motifs")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get reference sequence
    human_rDNA_ref = Path("references/human_rDNA.fa")
    seq = ""
    with open(human_rDNA_ref, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq += line.strip()
    
    # Extract the conserved ~400bp motif (first part of each unit)
    motif_seq = seq[13493:13493+400]
    
    # Create reference file
    motif_file = output_dir / "r_repeat_motif.fa"
    with open(motif_file, 'w') as f:
        f.write(f">r_repeat_motif\n{motif_seq}\n")
    
    logger.info(f"Created r-repeat motif reference in {motif_file}")
    return str(motif_file)

def align_to_motif(fastq_path, output_dir=None, reuse_existing=False):
    """
    Align r-repeat sequences to the conserved motif to identify units.
    
    Parameters:
        fastq_path: Path to FASTQ file created by samtools fastq
        output_dir: Directory for output files
        reuse_existing: Whether to reuse existing alignment if available
        
    Returns:
        str: Path to the alignment BAM file
    """
    # Ensure string paths
    fastq_path = str(fastq_path)
    
    if output_dir is None:
        output_dir = Path("temp_files")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get motif reference
    motif_ref = create_motif_reference()
    
    # Output BAM file
    alignment_path = output_dir / f"{Path(fastq_path).stem}_motif_aligned.bam"
    
    # Check if we can reuse existing alignment
    if reuse_existing and alignment_path.exists() and Path(f"{alignment_path}.bai").exists():
        logger.info(f"Using existing motif alignment: {alignment_path}")
        return str(alignment_path)
    
    # Run minimap2 with parameters optimized for detecting multiple motif occurrences
    # Convert all Path objects to strings
    cmd = [
        "minimap2",
        "-x", "map-ont",      # Oxford Nanopore preset
        "-k", "10",           # Shorter k-mer for higher sensitivity
        "-w", "5",            # Reduced window size
        "-a",                 # Output SAM format
        "-N", "100",          # Report up to 100 alignments per read
        "--secondary=yes",    # Report secondary alignments
        "-p", "0.5",          # Lower identity threshold for chaining (default 0.8)
        "-Y",                 # Use soft clipping for supplementary alignments
        str(motif_ref),       # Reference sequence - ensure string
        str(fastq_path)       # Query sequences (r-repeats) - ensure string
    ]
    
    logger.info(f"Aligning r-repeat sequences to motif reference with command: {' '.join(cmd)}")
    
    # Sort and convert to BAM
    with open(alignment_path, 'wb') as bam_file:
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ["samtools", "sort", "-o", str(alignment_path)],
            stdin=p1.stdout
        )
        p1.stdout.close()
        p2.communicate()
    
    # Index BAM file
    pysam.index(str(alignment_path))
    
    logger.info(f"Aligned r-repeat sequences to motif reference: {alignment_path}")
    return str(alignment_path)

def process_motif_alignments(alignment_bam):
    """
    Process minimap2 alignments to identify r-repeat units.
    
    This function reads the alignment file and identifies all instances
    of the conserved motif within each r-repeat region.
    
    Parameters:
        alignment_bam: Path to BAM file with motif alignments
        
    Returns:
        dict: Dictionary with r-repeat unit analysis results
    """
    
    # Group alignments by read/region
    read_alignments = {}
    
    # First pass: collect all alignments by read/region
    with pysam.AlignmentFile(alignment_bam, "rb") as bam:
        for aln in bam:
            if aln.is_unmapped:
                continue
            
            # Skip very low quality alignments immediately
            if aln.mapping_quality < 5:
                continue
                
            read_id = aln.query_name

            if read_id not in read_alignments:
                read_alignments[read_id] = []
            
            motif_alignment_len = aln.query_alignment_end - aln.query_alignment_start
            if motif_alignment_len < 250 or motif_alignment_len > 500:
                continue

            motif_length = 400
            if aln.has_tag('NM'):  # Number of mismatches
                nm = aln.get_tag('NM')
                identity = 1 - (nm / motif_length) if motif_length > 0 else 0
            else:
                identity = 0
            if identity < 0.5:
                continue
            
            # Store alignment - we now work with the reference orientation
            read_alignments[read_id].append(aln)
    
    logger.info(f"Processing motif alignments for {len(read_alignments)} r-repeat regions")
    
    return read_alignments

def correspond_units_to_regions(read_alignments, read_id_to_r_repeat):
    """
    Correspond r-repeat units to the original r-repeat regions.
    
    Args:
        read_alignments: Dictionary with read ID to motif alignments
        read_id_to_r_repeat: Dictionary with read ID to r-repeat regions created in r_repeat_detection.py
    """
    # Track unit count distribution
    unit_count_distribution = {}
    # Track unit lengths
    all_unit_lengths = []
    # Track r-repeat region lengths
    r_repeat_region_lengths = []
    # Keep track of read IDs with problematic units to exclude from final analysis
    read_ids_with_problematic_units = set()
    
    # Store comprehensive unit data for serialization
    unit_analysis_data = {
        "unit_counts": {},
        "units_by_read": {},
        "unit_length_stats": {},
        "problematic_units": {
            "short_units": [],
            "long_units": []
        }
    }
    
    for read_id, r_repeat_regions in read_id_to_r_repeat.items():
        # Handle both single region and multiple regions per read
        if isinstance(r_repeat_regions, list):
            regions_list = r_repeat_regions
        else:
            regions_list = [r_repeat_regions]
            
        for r_repeat_region in regions_list:
            # Add region length to the list for distribution analysis
            r_repeat_region_lengths.append(r_repeat_region['length'])
            
            region_units = []
            if read_id not in read_alignments:
                continue
            for alignment in read_alignments[read_id]:
                if alignment.query_alignment_start >= r_repeat_region['start'] and alignment.query_alignment_end <= r_repeat_region['end']:
                    region_units.append(alignment)
            if region_units:
                # Sort units by start position
                region_units = sorted(region_units, key=lambda a: a.query_alignment_start)
                r_repeat_region['units'] = region_units
                # Unit lengths are calculated by the distance between the start of the motif alignments
                unit_lengths = []
                has_problematic_unit = False
                
                if len(region_units) > 1:
                    for i in range(1, len(region_units)):
                        unit_length = region_units[i].query_alignment_start - region_units[i-1].query_alignment_start
                        unit_lengths.append(unit_length)
                        all_unit_lengths.append(unit_length)
                        
                        # Check if this unit is too long (>1000bp) or too short (<300bp)
                        if unit_length > 1000 or unit_length < 300:
                            has_problematic_unit = True
                            read_ids_with_problematic_units.add(read_id)
                
                # Add the last unit length (using distance from motif start to region end)
                last_unit_length = r_repeat_region['end'] - region_units[-1].query_alignment_start
                unit_lengths.append(last_unit_length)
                all_unit_lengths.append(last_unit_length)
                
                # Check if the last unit is too long or too short
                if last_unit_length > 1000 or last_unit_length < 300:
                    has_problematic_unit = True
                    read_ids_with_problematic_units.add(read_id)
                
                # Mark regions with problematic units for exclusion from final analysis
                r_repeat_region['has_problematic_unit'] = has_problematic_unit
                r_repeat_region['has_long_unit'] = any(length > 1000 for length in unit_lengths)
                r_repeat_region['has_short_unit'] = any(length < 300 for length in unit_lengths)
                
                r_repeat_region['unit_lengths'] = unit_lengths
                r_repeat_region['unit_count'] = len(region_units)
                
                # Only include regions without problematic units in the distribution calculations
                if not has_problematic_unit:
                    # Track unit count distribution
                    count = len(region_units)
                    if count not in unit_count_distribution:
                        unit_count_distribution[count] = 0
                    unit_count_distribution[count] += 1
        
        # Add to unit analysis data for serialization
        if isinstance(r_repeat_regions, list):
            if read_id not in unit_analysis_data["units_by_read"]:
                unit_analysis_data["units_by_read"][read_id] = []
            
            for i, region in enumerate(regions_list):
                if 'unit_count' in region and 'unit_lengths' in region:
                    region_data = {
                        "region_index": i,
                        "unit_count": region['unit_count'],
                        "unit_lengths": region['unit_lengths'],
                        "has_problematic_unit": region.get('has_problematic_unit', False),
                        "has_short_unit": region.get('has_short_unit', False),
                        "has_long_unit": region.get('has_long_unit', False)
                    }
                    unit_analysis_data["units_by_read"][read_id].append(region_data)
                    
                    # Track problematic units
                    if region.get('has_short_unit', False):
                        unit_analysis_data["problematic_units"]["short_units"].append({
                            "read_id": read_id,
                            "region_index": i,
                            "unit_lengths": region['unit_lengths']
                        })
                    if region.get('has_long_unit', False):
                        unit_analysis_data["problematic_units"]["long_units"].append({
                            "read_id": read_id,
                            "region_index": i,
                            "unit_lengths": region['unit_lengths']
                        })
        else:
            region = r_repeat_regions
            if 'unit_count' in region and 'unit_lengths' in region:
                region_data = {
                    "unit_count": region['unit_count'],
                    "unit_lengths": region['unit_lengths'],
                    "has_problematic_unit": region.get('has_problematic_unit', False),
                    "has_short_unit": region.get('has_short_unit', False),
                    "has_long_unit": region.get('has_long_unit', False)
                }
                unit_analysis_data["units_by_read"][read_id] = region_data
                
                # Track problematic units
                if region.get('has_short_unit', False):
                    unit_analysis_data["problematic_units"]["short_units"].append({
                        "read_id": read_id,
                        "unit_lengths": region['unit_lengths']
                    })
                if region.get('has_long_unit', False):
                    unit_analysis_data["problematic_units"]["long_units"].append({
                        "read_id": read_id,
                        "unit_lengths": region['unit_lengths']
                    })
    
    # Filter out regions with problematic units from the lengths lists
    filtered_unit_lengths = [length for length in all_unit_lengths if 300 <= length <= 1000]
    
    # Log unit count distribution
    total_regions = sum(unit_count_distribution.values())
    if total_regions > 0:
        logger.info(f"r-repeat unit count distribution across {total_regions} regions (excluding regions with units <300bp or >1000bp):")
        for count, num_regions in sorted(unit_count_distribution.items()):
            percentage = (num_regions / total_regions) * 100
            logger.info(f"  {count} units: {num_regions} regions ({percentage:.1f}%)")
    
    # Log unit length statistics
    if filtered_unit_lengths:
        logger.info(f"Unit length statistics (excluding problematic units): min={min(filtered_unit_lengths)}, median={np.median(filtered_unit_lengths):.1f}, "
                   f"mean={np.mean(filtered_unit_lengths):.1f}, max={max(filtered_unit_lengths)}")
    
    # Count and log excluded reads with short and long units
    short_unit_reads = sum(1 for read_id, regions in read_id_to_r_repeat.items() 
                          if isinstance(regions, list) 
                          and any(region.get('has_short_unit', False) for region in regions)
                          or (not isinstance(regions, list) and regions.get('has_short_unit', False)))
                          
    long_unit_reads = sum(1 for read_id, regions in read_id_to_r_repeat.items() 
                          if isinstance(regions, list) 
                          and any(region.get('has_long_unit', False) for region in regions)
                          or (not isinstance(regions, list) and regions.get('has_long_unit', False)))
    
    logger.info(f"Found and marked {len(read_ids_with_problematic_units)} reads with problematic units for exclusion from analysis:")
    logger.info(f"  - {short_unit_reads} reads with units <300bp")
    logger.info(f"  - {long_unit_reads} reads with units >1000bp")
    
    # Log region length statistics
    if r_repeat_region_lengths:
        logger.info(f"r-repeat region length statistics: min={min(r_repeat_region_lengths)}, median={np.median(r_repeat_region_lengths):.1f}, "
                   f"mean={np.mean(r_repeat_region_lengths):.1f}, max={max(r_repeat_region_lengths)}")
    
    # Update unit count distribution in serializable format
    for count, freq in unit_count_distribution.items():
        unit_analysis_data["unit_counts"][str(count)] = freq
    
    # Add statistics for unit lengths
    if filtered_unit_lengths:
        unit_analysis_data["unit_length_stats"] = {
            "count": len(filtered_unit_lengths),
            "min": min(filtered_unit_lengths),
            "max": max(filtered_unit_lengths),
            "mean": float(np.mean(filtered_unit_lengths)),
            "median": float(np.median(filtered_unit_lengths)),
            "std_dev": float(np.std(filtered_unit_lengths)),
            "quartiles": {
                "q1": float(np.percentile(filtered_unit_lengths, 25)),
                "q2": float(np.percentile(filtered_unit_lengths, 50)),  # Same as median
                "q3": float(np.percentile(filtered_unit_lengths, 75))
            }
        }
    
    # Return both the filtered and unfiltered data, plus read IDs to exclude
    return unit_count_distribution, filtered_unit_lengths, r_repeat_region_lengths, read_ids_with_problematic_units, unit_analysis_data

def generate_unit_count_visualization(unit_count_distribution, output_dir=None):
    """
    Generate visualization of unit count distribution.
    
    Parameters:
        unit_count_distribution: Dictionary mapping unit counts to number of regions
        output_dir: Directory for output files
        
    Returns:
        str: Path to visualization file
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if not unit_count_distribution:
        logger.warning("No unit count data to visualize")
        return None
        
    if output_dir is None:
        output_dir = Path("results/r_repeat_units")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Prepare data
    counts = []
    freqs = []
    for count, freq in sorted(unit_count_distribution.items()):
        counts.append(count)
        freqs.append(freq)
    
    # Create visualization
    plt.figure(figsize=(10, 6))
    
    # Bar plot with counts
    bars = plt.bar(counts, freqs, alpha=0.7)
    
    # Add count labels on top of bars
    for bar, freq in zip(bars, freqs):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{freq}', ha='center', va='bottom')
    
    plt.title("Distribution of r-repeat Units per Region", fontsize=14)
    plt.xlabel("Number of Units", fontsize=12)
    plt.ylabel("Number of r-repeat Regions", fontsize=12)
    plt.xticks(counts)
    plt.grid(axis='y', alpha=0.3)
    
    # Save figure
    output_file = output_dir / "unit_count_distribution.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"Saved unit count distribution visualization to {output_file}")
    return str(output_file)

def generate_unit_length_visualization(unit_lengths, output_dir=None):
    """
    Generate visualization of unit length distribution.
    
    Parameters:
        unit_lengths: List of unit lengths
        output_dir: Directory for output files
        
    Returns:
        str: Path to visualization file
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if not unit_lengths:
        logger.warning("No unit length data to visualize")
        return None
        
    if output_dir is None:
        output_dir = Path("results/r_repeat_units")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create visualization
    plt.figure(figsize=(10, 6))
    
    # Histogram plot with density curve
    sns.histplot(unit_lengths, bins=30, kde=True)
    
    # Add reference line for expected ~600-700bp unit length
    plt.axvline(x=650, color='r', linestyle='--', label='Expected typical unit length (~650bp)')
    
    plt.title("Distribution of r-repeat Unit Lengths", fontsize=14)
    plt.xlabel("Unit Length (bp)", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    
    # Add statistics in text box
    stats_text = f"n = {len(unit_lengths)}\n"
    stats_text += f"min = {min(unit_lengths):.0f} bp\n"
    stats_text += f"mean = {np.mean(unit_lengths):.1f} bp\n"
    stats_text += f"median = {np.median(unit_lengths):.1f} bp\n"
    stats_text += f"max = {max(unit_lengths):.0f} bp"
    
    plt.annotate(stats_text, xy=(0.02, 0.95), xycoords='axes fraction',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
                 va='top', fontsize=10)
    
    # Save figure
    output_file = output_dir / "unit_length_distribution.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"Saved unit length distribution visualization to {output_file}")
    return str(output_file)

def generate_region_length_visualization(region_lengths, output_dir=None):
    """
    Generate visualization of r-repeat region length distribution.
    
    Parameters:
        region_lengths: List of r-repeat region lengths
        output_dir: Directory for output files
        
    Returns:
        str: Path to visualization file
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if not region_lengths:
        logger.warning("No region length data to visualize")
        return None
        
    if output_dir is None:
        output_dir = Path("results/r_repeat_units")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create visualization
    plt.figure(figsize=(10, 6))
    
    # Histogram plot with density curve
    sns.histplot(region_lengths, bins=30, kde=True)
    
    plt.title("Distribution of r-repeat Region Lengths", fontsize=14)
    plt.xlabel("Region Length (bp)", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.grid(axis='y', alpha=0.3)
    
    # Add statistics in text box
    stats_text = f"n = {len(region_lengths)}\n"
    stats_text += f"min = {min(region_lengths):.0f} bp\n"  # Removed extra colon
    stats_text += f"mean = {np.mean(region_lengths):.1f} bp\n"
    stats_text += f"median = {np.median(region_lengths):.1f} bp\n"  # Fixed extra colon in format specifier
    stats_text += f"max = {max(region_lengths):.0f} bp"
    
    plt.annotate(stats_text, xy=(0.02, 0.95), xycoords='axes fraction',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
                 va='top', fontsize=10)
    
    # Save figure
    output_file = output_dir / "region_length_distribution.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"Saved r-repeat region length distribution visualization to {output_file}")
    return str(output_file)

def find_long_units(read_alignments, read_id_to_r_repeat, threshold=1000):
    """
    Find r-repeat regions with unusually long units (>threshold).
    
    Parameters:
        read_alignments: Dictionary with read ID to motif alignments
        read_id_to_r_repeat: Dictionary with read ID to r-repeat regions
        threshold: Length threshold in bp (default: 1000)
        
    Returns:
        list: List of dictionaries containing read_id, region_index, and unit length info
    """
    long_units = []
    
    for read_id, r_repeat_regions in read_id_to_r_repeat.items():
        if read_id not in read_alignments:
            continue
            
        # Handle both single region and multiple regions per read
        if isinstance(r_repeat_regions, list):
            regions_list = r_repeat_regions
        else: regions_list = [r_repeat_regions]
        
        for region_idx, r_repeat_region in enumerate(regions_list):
            region_units = []
            
            # Find alignments for this region
            for alignment in read_alignments[read_id]:
                if (alignment.query_alignment_start >= r_repeat_region['start'] and 
                    alignment.query_alignment_end <= r_repeat_region['end']):
                    region_units.append(alignment)
            
            if len(region_units) > 1:
                # Sort units by start position
                region_units = sorted(region_units, key=lambda a: a.query_alignment_start)
                
                # Calculate unit lengths
                for i in range(1, len(region_units)):
                    unit_length = region_units[i].query_alignment_start - region_units[i-1].query_alignment_start
                    
                    # Check if unit exceeds threshold
                    if unit_length > threshold:
                        # Add to results
                        long_units.append({
                            'read_id': read_id,
                            'region_index': region_idx,
                            'unit_index': i-1,  # 0-based index of unit
                            'unit_length': unit_length,
                            'start': region_units[i-1].query_alignment_start,
                            'end': region_units[i].query_alignment_start
                        })
                
                # Check last unit
                if len(region_units) > 0:
                    last_unit_length = r_repeat_region['end'] - region_units[-1].query_alignment_start
                    if last_unit_length > threshold:
                        long_units.append({
                            'read_id': read_id,
                            'region_index': region_idx,
                            'unit_index': len(region_units)-1,  # 0-based index of unit
                            'unit_length': last_unit_length,
                            'start': region_units[-1].query_alignment_start,
                            'end': r_repeat_region['end']
                        })
    
    # Sort by unit length (descending)
    long_units.sort(key=lambda x: x['unit_length'], reverse=True)
    
    logger.info(f"Found {len(long_units)} units longer than {threshold}bp")
    
    # Print the results to log
    if long_units:
        logger.info("Top 10 longest units:")
        for i, unit in enumerate(long_units[:10]):
            logger.info(f"{i+1}. Read {unit['read_id']} (region {unit['region_index']}, "
                      f"unit {unit['unit_index']}): {unit['unit_length']}bp")
    
    print(long_units)
    return long_units

def generate_unit_visualizations(unit_analysis_data, output_dir=None):
    """
    Generate all unit visualizations and save serialized data for cross-sample comparison.
    
    Parameters:
        unit_analysis_data: Dictionary with comprehensive unit analysis data
        output_dir: Directory for output files
        
    Returns:
        dict: Dictionary mapping visualization types to file paths
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if output_dir is None:
        output_dir = Path("r_repeat_results/unit_analysis")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create stats directory for serialized visualization data
    stats_dir = output_dir / "stats"
    stats_dir.mkdir(exist_ok=True)
    
    # Dictionary to store paths to visualizations
    visualizations = {}
    
    # 1. Unit count distribution visualization
    if "unit_counts" in unit_analysis_data and unit_analysis_data["unit_counts"]:
        unit_counts = {int(k): v for k, v in unit_analysis_data["unit_counts"].items()}
        viz_path = generate_unit_count_visualization(unit_counts, output_dir)
        if viz_path:
            visualizations["unit_count_distribution"] = viz_path
    
    # 2. Unit length distribution visualization
    if ("unit_length_stats" in unit_analysis_data and 
        "filtered_unit_lengths" in unit_analysis_data and 
        unit_analysis_data["filtered_unit_lengths"]):
        viz_path = generate_unit_length_visualization(
            unit_analysis_data["filtered_unit_lengths"], output_dir)
        if viz_path:
            visualizations["unit_length_distribution"] = viz_path
    
    # 3. Region length distribution visualization
    if "region_lengths" in unit_analysis_data and unit_analysis_data["region_lengths"]:
        viz_path = generate_region_length_visualization(
            unit_analysis_data["region_lengths"], output_dir)
        if viz_path:
            visualizations["region_length_distribution"] = viz_path
    
    # Save serialized statistics for cross-sample comparison
    with open(stats_dir / "unit_analysis_stats.json", 'w') as f:
        # Extract the key statistics for cross-sample comparison
        comparison_stats = {
            "unit_counts": unit_analysis_data.get("unit_counts", {}),
            "unit_length_stats": unit_analysis_data.get("unit_length_stats", {}),
            "problematic_units": {
                "short_unit_count": len(unit_analysis_data.get("problematic_units", {}).get("short_units", [])),
                "long_unit_count": len(unit_analysis_data.get("problematic_units", {}).get("long_units", []))
            }
        }
        json.dump(comparison_stats, f, indent=2)
    
    logger.info(f"Generated {len(visualizations)} unit visualizations in {output_dir}")
    logger.info(f"Saved serialized unit statistics to {stats_dir}")
    
    return visualizations

# Sample usage (to be expanded with actual unit identification logic)
if __name__ == "__main__":
    # This is just a test to verify the basic functions work
    logger.info("Testing r-repeat unit analysis functions")
    
    # Create reference
    motif_ref = create_motif_reference()
    logger.info(f"Created motif reference at {motif_ref}")
    
    # Test reverse complement function
    test_seq = "ACGTACGT"
    rev_comp = reverse_complement(test_seq)
    logger.info(f"Reverse complement of {test_seq} is {rev_comp}")
