"""
Methylation analysis module for r-repeat regions.

This module focuses on extracting and analyzing methylation patterns in 
previously identified r-repeat regions, with special attention to unit structure.
It also includes analysis of upstream coding region methylation patterns.
"""

import pickle
import logging
import pysam
import numpy as np
import sys
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Direct imports from r_repeat package
from r_repeat.upstream_analysis.analyze_upstream_coding import analyze_upstream_coding_methylation
from r_repeat.downstream_analysis.analyze_downstream_igs import analyze_downstream_igs_methylation

# Ensure parent directory is in path for imports
parent_dir = str(Path(__file__).parent.parent.parent)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

# Try importing modules from reorganized structure first
try:
    from ..upstream_analysis.analyze_upstream_coding import analyze_upstream_coding_methylation
    from ..downstream_analysis.analyze_downstream_igs import analyze_downstream_igs_methylation
except ImportError:
    # Fallback to original import paths
    try:
        from analyze_upstream_coding import analyze_upstream_coding_methylation
        try:
            from analyze_downstream_igs import analyze_downstream_igs_methylation
        except ImportError:
            analyze_downstream_igs_methylation = None
            print("Warning: analyze_downstream_igs module not found. Downstream analysis will be disabled.")
    except ImportError:
        analyze_upstream_coding_methylation = None
        analyze_downstream_igs_methylation = None
        print("Warning: analyze_upstream_coding module not found. Upstream and downstream analysis will be disabled.")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler("logs/r_repeat_methylation.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define a fallback function if analyze_downstream_igs is not available
def fallback_downstream_analysis(*args, **kwargs):
    """
    Fallback function when analyze_downstream_igs module is not available.
    
    Returns:
        dict: Minimal result structure with error status
    """
    logger.warning("analyze_downstream_igs module not found. Downstream analysis skipped.")
    return {
        'status': 'error',
        'message': 'analyze_downstream_igs module not found. Please ensure the file exists and is importable.',
        'visualizations': {}
    }

def extract_r_repeat_methylation(r_repeat_regions, mod_bam_path, output_dir=None):
    """
    Extract methylation data for identified r-repeat regions.
    
    Parameters:
        r_repeat_regions: Dictionary mapping read IDs to r-repeat regions (single or multiple per read)
        mod_bam_path: Path to BAM file with modification (methylation) data
        output_dir: Directory for output files
        
    Returns:
        dict: Dictionary with r-repeat methylation data
    """
    # Convert input paths to Path objects
    mod_bam_path = Path(mod_bam_path)
    sample_name = mod_bam_path.stem
    
    if output_dir is None:
        output_dir = Path("r_repeat_results") / sample_name / "methylation"
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create data directory for serialized data
    data_dir = output_dir.parent / "data" 
    data_dir.mkdir(exist_ok=True, parents=True)
    
    # Dictionary to store methylation data
    r_repeat_data = {}
    
    # Check if we have multiple regions per read
    has_multiple_regions = any(isinstance(regions, list) for regions in r_repeat_regions.values())
    
    # Extract methylation data from BAM
    with pysam.AlignmentFile(str(mod_bam_path), "rb", check_sq=False) as bam:
        for read in bam:
            read_id = read.query_name
            if read_id not in r_repeat_regions:
                continue
            
            # Handle either single or multiple regions
            if has_multiple_regions:
                # Multiple regions per read
                r_repeat_data[read_id] = []
                regions = r_repeat_regions[read_id]
                
                for region_idx, coords in enumerate(regions):
                    r_start = coords['start']
                    r_end = coords['end']
                    r_length = coords['length']
                    is_reverse = coords['is_reverse']
                    
                    # Extract methylation data for this region
                    methylation_data = extract_methylation_for_region(read, r_start, r_end, r_length, is_reverse)
                    
                    # Extract methylation data for coding region before this r-repeat region
                    upstream_methylation = extract_upstream_coding_methylation(read, r_start, is_reverse)
                    
                    # Extract methylation data for IGS region after this r-repeat region
                    downstream_methylation = extract_downstream_igs_methylation(read, r_end, is_reverse)
                    
                    # Get unit information if available
                    unit_data = None
                    if 'units' in coords:
                        unit_data = extract_unit_methylation(coords['units'], methylation_data)
                        # Remove the units from coords to avoid pickling errors
                        coords_copy = coords.copy()
                        if 'units' in coords_copy:
                            del coords_copy['units']
                    else:
                        coords_copy = coords
                    
                    if methylation_data:
                        region_data = {
                            'coordinates': coords_copy,
                            'methylation_data': sorted(methylation_data),
                            'avg_methylation': np.mean([score for _, score in methylation_data]),
                            'methylation_count': len(methylation_data),
                            'methylation_density': len(methylation_data) / r_length if r_length > 0 else 0,
                            'upstream_coding_methylation': upstream_methylation,
                            'downstream_igs_methylation': downstream_methylation,
                            'region_index': region_idx,
                            'unit_methylation': unit_data
                        }
                        r_repeat_data[read_id].append(region_data)
            else:
                # Single region per read (original behavior)
                coords = r_repeat_regions[read_id]
                r_start = coords['start']
                r_end = coords['end']
                r_length = coords['length']
                is_reverse = coords['is_reverse']
                
                # Extract methylation data for this region
                methylation_data = extract_methylation_for_region(read, r_start, r_end, r_length, is_reverse)
                
                # Extract methylation data for coding region before this r-repeat region
                upstream_methylation = extract_upstream_coding_methylation(read, r_start, is_reverse)
                
                # Extract methylation data for IGS region after this r-repeat region
                downstream_methylation = extract_downstream_igs_methylation(read, r_end, is_reverse)
                
                # Get unit information if available
                unit_data = None
                if 'units' in coords:
                    unit_data = extract_unit_methylation(coords['units'], methylation_data)
                    # Remove the units from coords to avoid pickling errors
                    coords_copy = coords.copy()
                    if 'units' in coords_copy:
                        del coords_copy['units']
                else:
                    coords_copy = coords
                
                if methylation_data:
                    r_repeat_data[read_id] = {
                        'coordinates': coords_copy,
                        'methylation_data': sorted(methylation_data),
                        'avg_methylation': np.mean([score for _, score in methylation_data]),
                        'methylation_count': len(methylation_data),
                        'methylation_density': len(methylation_data) / r_length if r_length > 0 else 0,
                        'upstream_coding_methylation': upstream_methylation,
                        'downstream_igs_methylation': downstream_methylation,
                        'unit_methylation': unit_data
                    }
    
    # Save methylation data to both directories for consistency
    output_file = data_dir / "r_repeat_methylation.pkl"
    with open(output_file, "wb") as f:
        pickle.dump(r_repeat_data, f)
    
    # Also save a metadata JSON with basic stats for easier cross-sample comparison
    metadata = {
        "sample_name": sample_name,
        "extraction_date": str(pd.Timestamp.now()),
        "region_count": (
            sum(len(regions) for regions in r_repeat_data.values()) 
            if has_multiple_regions else len(r_repeat_data)
        ),
        "read_count": len(r_repeat_data),
        "average_methylation": float(
            np.mean([
                reg['avg_methylation'] for regions in r_repeat_data.values()
                for reg in (regions if has_multiple_regions else [regions])
            ]) / 255
        ) if r_repeat_data else 0.0
    }
    
    # Save metadata
    with open(data_dir / "methylation_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    # Log different statistics based on structure
    if has_multiple_regions:
        region_count = sum(len(regions) for regions in r_repeat_data.values())
        logger.info(f"Extracted methylation data for {region_count} r-repeat regions from {len(r_repeat_data)} reads")
    else:
        logger.info(f"Extracted methylation data for {len(r_repeat_data)} r-repeat regions")
    
    return r_repeat_data

def extract_methylation_for_region(read, r_start, r_end, r_length, is_reverse):
    """
    Helper function to extract methylation data for a specific region
    
    Parameters:
        read: pysam alignment record with methylation data
        r_start: start position of r-repeat region
        r_end: end position of r-repeat region
        r_length: length of r-repeat region
        is_reverse: whether the region is on the reverse strand
        
    Returns:
        list: List of (position, methylation_score) tuples
    """
    methylation_data = []
    if hasattr(read, 'modified_bases') and read.modified_bases:
        # Find the 5mC methylation key (typically ('C', 0, 'm') for CpG methylation)
        meth_key = None
        for key in read.modified_bases.keys():
            if key[0] == 'C' and key[2] == 'm':
                meth_key = key
                break
        
        if meth_key:
            # Extract methylation positions and scores within the r-repeat
            base_mods = read.modified_bases[meth_key]
            
            # Handle strand orientation - flip methylation positions for reverse strand 
            # instead of flipping the sequence
            if is_reverse:
                # For reverse strand, we need to flip the methylation positions
                read_len = read.query_length
                base_mods = [(read_len - pos - 2, score) for pos, score in base_mods]
            
            # Filter for positions within the r-repeat region
            for pos, score in base_mods:
                if r_start <= pos <= r_end:
                    # Normalize position to be relative to r-repeat start
                    norm_pos = pos - r_start
                    methylation_data.append((norm_pos, score))
    
    return methylation_data

def extract_upstream_coding_methylation(read, r_start, is_reverse):
    """
    Extract methylation data for the coding region upstream of an r-repeat region.
    
    Parameters:
        read: pysam alignment record with methylation data
        r_start: start position of r-repeat region
        is_reverse: whether the region is on the reverse strand
        
    Returns:
        dict: Dictionary with upstream coding region methylation information
    """
    # Define the upstream region (2000bp before r-repeat start)
    upstream_length = 2000
    
    # Check if there's enough space for upstream region
    if r_start < upstream_length:
        # Not enough space, return -1 to signal unavailability
        return {
            'avg_methylation': -1,
            'methylation_count': 0,
            'methylation_data': [],
            'available': False
        }
    
    upstream_start = r_start - upstream_length
    upstream_end = r_start
    
    # Extract methylation data
    methylation_data = []
    if hasattr(read, 'modified_bases') and read.modified_bases:
        # Find the 5mC methylation key (typically ('C', 0, 'm') for CpG methylation)
        meth_key = None
        for key in read.modified_bases.keys():
            if key[0] == 'C' and key[2] == 'm':
                meth_key = key
                break
        
        if meth_key:
            # Extract methylation positions and scores within the upstream region
            base_mods = read.modified_bases[meth_key]
            
            # Handle strand orientation - flip methylation positions for reverse strand
            if is_reverse:
                read_len = read.query_length
                base_mods = [(read_len - pos - 2, score) for pos, score in base_mods]
            
            # Filter for positions within the upstream coding region
            for pos, score in base_mods:
                if upstream_start <= pos < upstream_end:
                    # Normalize position to be relative to upstream region start
                    norm_pos = pos - upstream_start
                    methylation_data.append((norm_pos, score))
    
    if methylation_data:
        avg_methylation = np.mean([score for _, score in methylation_data])
        return {
            'avg_methylation': avg_methylation,
            'methylation_count': len(methylation_data),
            'methylation_data': sorted(methylation_data),
            'methylation_density': len(methylation_data) / upstream_length,
            'available': True
        }
    else:
        # No methylation data found in upstream region
        return {
            'avg_methylation': 0,
            'methylation_count': 0,
            'methylation_data': [],
            'methylation_density': 0,
            'available': True
        }

def extract_downstream_igs_methylation(read, r_end, is_reverse):
    """
    Extract methylation data for the intergenic spacer (IGS) region downstream of an r-repeat region.
    
    Parameters:
        read: pysam alignment record with methylation data
        r_end: end position of r-repeat region
        is_reverse: whether the region is on the reverse strand
        
    Returns:
        dict: Dictionary with downstream IGS region methylation information
    """
    # Define the downstream region (2000bp after r-repeat end)
    downstream_length = 2000
    
    # Define the downstream region boundaries
    downstream_start = r_end
    downstream_end = r_end + downstream_length
    
    # Extract methylation data
    methylation_data = []
    if hasattr(read, 'modified_bases') and read.modified_bases:
        # Find the 5mC methylation key (typically ('C', 0, 'm') for CpG methylation)
        meth_key = None
        for key in read.modified_bases.keys():
            if key[0] == 'C' and key[2] == 'm':
                meth_key = key
                break
        
        if meth_key:
            # Extract methylation positions and scores within the downstream region
            base_mods = read.modified_bases[meth_key]
            
            # Handle strand orientation - flip methylation positions for reverse strand
            if is_reverse:
                read_len = read.query_length
                base_mods = [(read_len - pos - 2, score) for pos, score in base_mods]
            
            # Filter for positions within the downstream IGS region
            for pos, score in base_mods:
                if downstream_start <= pos < downstream_end:
                    # Normalize position to be relative to downstream region start
                    norm_pos = pos - downstream_start
                    methylation_data.append((norm_pos, score))
    
    if methylation_data:
        avg_methylation = np.mean([score for _, score in methylation_data])
        return {
            'avg_methylation': avg_methylation,
            'methylation_count': len(methylation_data),
            'methylation_data': sorted(methylation_data),
            'methylation_density': len(methylation_data) / downstream_length,
            'available': True
        }
    else:
        # No methylation data found in downstream region
        return {
            'avg_methylation': 0,
            'methylation_count': 0,
            'methylation_data': [],
            'methylation_density': 0,
            'available': len(read.query_sequence) >= r_end + downstream_length  # Only mark as available if region exists
        }

def extract_unit_methylation(units, methylation_data):
    """
    Extract methylation data for each unit in an r-repeat region.
    
    Parameters:
        units: List of unit objects with alignment information
        methylation_data: List of (position, methylation_score) tuples
        
    Returns:
        dict: Dictionary with unit-specific methylation data
    """
    unit_methylation = {}
    
    # Skip if no units or methylation data
    if not units or not methylation_data:
        return unit_methylation
        
    # Sort units by position - extract the query_alignment_start value directly to avoid pickling issues
    sorted_units = sorted([(i, u.query_alignment_start, u.query_alignment_end) 
                          for i, u in enumerate(units)], 
                          key=lambda x: x[1])
    
    # Process each unit
    for i in range(len(sorted_units)):
        unit_idx, unit_start, unit_end = sorted_units[i]
        # Determine the end boundary - either the next unit's start or the current unit's end
        next_unit_start = sorted_units[i+1][1] if i < len(sorted_units) - 1 else unit_end
        
        # Find methylation data for this unit
        unit_meth = []
        for pos, score in methylation_data:
            if unit_start <= pos < next_unit_start:
                # Normalize position relative to unit start
                norm_pos = pos - unit_start
                unit_meth.append((norm_pos, score))
        
        if unit_meth:
            unit_methylation[f"unit_{unit_idx+1}"] = {
                'methylation_data': unit_meth,
                'avg_methylation': np.mean([score for _, score in unit_meth]),
                'methylation_count': len(unit_meth),
                'unit_start': unit_start,
                'unit_end': next_unit_start,
                'unit_length': next_unit_start - unit_start
            }
    
    return unit_methylation

def visualize_methylation_patterns(r_repeat_data, output_dir=None):
    """
    Create visualizations of r-repeat methylation patterns with a focus on unit structure.
    
    Parameters:
        r_repeat_data: Dictionary with r-repeat methylation data (single or multiple regions per read)
        output_dir: Directory for output files
        
    Returns:
        list: List of paths to generated visualization files
    """
    # Set up visualization directory
    if output_dir is None:
        # Try to determine sample name from data
        sample_name = "unknown_sample"
        output_dir = Path("r_repeat_results") / sample_name / "methylation" / "figures"
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a stats directory for serialized visualization data
    stats_dir = output_dir.parent / "stats"
    stats_dir.mkdir(exist_ok=True)
    
    output_files = []
    visualization_stats = {}
    
    # Skip visualization if no data
    if not r_repeat_data:
        logger.warning("No methylation data to visualize")
        return output_files
    
    # Check if we have multiple regions per read
    has_multiple_regions = any(isinstance(regions, list) for regions in r_repeat_data.values())
    
    # Flatten data if needed to simplify visualization code, filtering out regions with problematic units
    flattened_data = []
    if has_multiple_regions:
        for read_id, regions in r_repeat_data.items():
            for region_data in regions:
                # Only include regions without problematic units in statistical visualizations
                if not region_data.get('has_problematic_unit', False):
                    flattened_data.append(region_data)
    else:
        # Include only regions without problematic units in visualizations
        flattened_data = [data for data in r_repeat_data.values() 
                         if not data.get('has_problematic_unit', False)]
    
    logger.info(f"Using {len(flattened_data)} regions for statistical visualization (excluding regions with problematic units)")
    
    # Figure 1: Methylation by unit position
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Collect unit methylation data
    unit_avg_meth = []
    for data in flattened_data:
        if 'unit_methylation' in data and data['unit_methylation']:
            for unit_id, unit_data in data['unit_methylation'].items():
                unit_num = int(unit_id.split('_')[1])  # Extract unit number from "unit_X"
                unit_avg_meth.append({
                    'unit_num': unit_num,
                    'avg_methylation': unit_data['avg_methylation'] / 255  # Normalize to 0-1
                })
    
    if unit_avg_meth:
        # Create DataFrame for easier plotting
        import pandas as pd
        df = pd.DataFrame(unit_avg_meth)
        
        # Box plot by unit position
        sns.boxplot(x='unit_num', y='avg_methylation', data=df, ax=ax)
        ax.set_title("Methylation by Unit Position", fontsize=14)
        ax.set_xlabel("Unit Number", fontsize=12)
        ax.set_ylabel("Average Methylation (0-1)", fontsize=12)
        ax.set_ylim(0, 1)
        
        unit_meth_file = output_dir / "methylation_by_unit.png"
        plt.savefig(unit_meth_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(unit_meth_file))
    
    # Figure 2: Average methylation by unit count
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Group by unit count
    unit_count_meth = {}
    for data in flattened_data:
        if 'coordinates' in data and 'unit_count' in data['coordinates']:
            unit_count = data['coordinates']['unit_count']
            avg_meth = data['avg_methylation'] / 255  # Normalize to 0-1
            
            if unit_count not in unit_count_meth:
                unit_count_meth[unit_count] = []
            unit_count_meth[unit_count].append(avg_meth)
    
    if unit_count_meth:
        # Prepare data for plotting
        x = []
        y = []
        for count in sorted(unit_count_meth.keys()):
            x.append(count)
            y.append(np.mean(unit_count_meth[count]))
        
        # Bar plot
        ax.bar(x, y, alpha=0.7)
        ax.set_title("Average Methylation by Unit Count", fontsize=14)
        ax.set_xlabel("Number of Units", fontsize=12)
        ax.set_ylabel("Average Methylation (0-1)", fontsize=12)
        ax.set_xticks(x)
        ax.set_ylim(0, 1)
        
        count_meth_file = output_dir / "methylation_by_unit_count.png"
        plt.savefig(count_meth_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(count_meth_file))
    
    # Figure 3: Averaged methylation patterns across r-repeat region length, classified by unit count
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Group r-repeats by unit count
    unit_counts = sorted(set(data['coordinates'].get('unit_count', 0) for data in flattened_data 
                          if 'coordinates' in data and 'unit_count' in data['coordinates'] 
                          and data['coordinates']['unit_count'] > 0))
    
    if unit_counts:
        # Create a line for each unit count
        for unit_count in unit_counts:
            # Get all r-repeats with this unit count
            regions = [data for data in flattened_data 
                      if 'coordinates' in data 
                      and 'unit_count' in data['coordinates'] 
                      and data['coordinates']['unit_count'] == unit_count]
            
            if not regions:
                continue
                
            # Group methylation data by relative position
            # Convert to relative position (0-100%) to normalize for different region lengths
            all_positions = {}
            for region in regions:
                length = region['coordinates']['length']
                for pos, score in region['methylation_data']:
                    rel_pos = int((pos / length) * 100)  # Convert to 0-100% position
                    if rel_pos not in all_positions:
                        all_positions[rel_pos] = []
                    all_positions[rel_pos].append(score / 255)  # Normalize score to 0-1
            
            # Calculate average methylation at each relative position
            positions = []
            avg_scores = []
            for pos in sorted(all_positions.keys()):
                positions.append(pos)
                avg_scores.append(np.mean(all_positions[pos]))
            
            # Smooth the plot with a rolling average
            window = 3  # Window size for rolling average
            if len(avg_scores) > window:
                smoothed_scores = np.convolve(avg_scores, np.ones(window)/window, mode='valid')
                # Adjust positions to match smoothed_scores length
                smoothed_positions = positions[window//2:-(window//2) if window//2 > 0 else None]
                
                ax.plot(smoothed_positions, smoothed_scores, 
                       label=f"{unit_count} units (n={len(regions)})",
                       linewidth=2, alpha=0.8)
            else:
                ax.plot(positions, avg_scores, 
                       label=f"{unit_count} units (n={len(regions)})",
                       linewidth=2, alpha=0.8)
        
        ax.set_title("Methylation Patterns Across r-repeat Regions by Unit Count", fontsize=14)
        ax.set_xlabel("Relative Position in r-repeat Region (%)", fontsize=12)
        ax.set_ylabel("Average Methylation (0-1)", fontsize=12)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend(title="Unit count", loc="best")
        
        pattern_file = output_dir / "methylation_pattern_by_unit_count.png"
        plt.savefig(pattern_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(pattern_file))
    
    # Figure 4: Methylation heatmap across units
    # Create a heatmap showing methylation patterns across individual units within r-repeats
    unit_methylation_data = {}
    for data in flattened_data:
        if 'unit_methylation' in data and data['unit_methylation']:
            for unit_id, unit_data in data['unit_methylation'].items():
                unit_num = int(unit_id.split('_')[1])  # Extract unit number
                if unit_num not in unit_methylation_data:
                    unit_methylation_data[unit_num] = []
                
                # Get average methylation for this unit
                unit_methylation_data[unit_num].append(unit_data['avg_methylation'] / 255)
    
    if unit_methylation_data and len(unit_methylation_data) > 1:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Convert to format suitable for heatmap
        import pandas as pd
        heatmap_data = []
        for unit_num in sorted(unit_methylation_data.keys()):
            values = unit_methylation_data[unit_num]
            for i, val in enumerate(values):
                heatmap_data.append({
                    'Unit Number': f"Unit {unit_num}",
                    'Sample': i,
                    'Methylation': val
                })
        
        df = pd.DataFrame(heatmap_data)
        
        # Create pivot table for heatmap
        pivot_df = df.pivot_table(
            index='Sample', 
            columns='Unit Number',
            values='Methylation',
            aggfunc='mean'
        )
        
        # Plot heatmap for units with sufficient data
        if not pivot_df.empty and pivot_df.shape[1] > 1:
            sns.heatmap(pivot_df, cmap="viridis", cbar_kws={'label': 'Methylation Level (0-1)'})
            plt.title("Methylation Patterns Across r-repeat Units", fontsize=14)
            
            heatmap_file = output_dir / "unit_methylation_heatmap.png"
            plt.savefig(heatmap_file, dpi=300, bbox_inches="tight")
            plt.close()
            output_files.append(str(heatmap_file))
    
    # Add stats for the visualizations
    if unit_avg_meth:
        visualization_stats["methylation_by_unit"] = {
            "unit_positions": sorted(set(item["unit_num"] for item in unit_avg_meth)),
            "methylation_by_unit": {
                unit: np.mean([item["avg_methylation"] for item in unit_avg_meth if item["unit_num"] == unit])
                for unit in set(item["unit_num"] for item in unit_avg_meth)
            },
            "sample_count": len(unit_avg_meth)
        }
    
    if unit_count_meth:
        visualization_stats["methylation_by_unit_count"] = {
            "unit_counts": list(unit_count_meth.keys()),
            "average_methylation": {
                count: np.mean(values) for count, values in unit_count_meth.items()
            },
            "sample_count": sum(len(values) for values in unit_count_meth.values())
        }
    
    # Add a new figure: Upstream coding methylation vs r-repeat methylation
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Collect data for the plot
    upstream_meth = []
    r_repeat_meth = []
    for data in flattened_data:
        # Check if upstream methylation data is available
        if 'upstream_coding_methylation' in data and data['upstream_coding_methylation']['available']:
            upstream_avg = data['upstream_coding_methylation']['avg_methylation'] / 255  # Normalize to 0-1
            r_repeat_avg = data['avg_methylation'] / 255  # Normalize to 0-1
            
            # Only include points where upstream methylation could be calculated
            if upstream_avg >= 0:
                upstream_meth.append(upstream_avg)
                r_repeat_meth.append(r_repeat_avg)
    
    if upstream_meth and r_repeat_meth:
        # Scatter plot with regression line
        ax.scatter(upstream_meth, r_repeat_meth, alpha=0.6)
        
        # Add regression line
        if len(upstream_meth) > 1:  # Need at least 2 points for regression
            from scipy.stats import linregress
            slope, intercept, r_value, p_value, std_err = linregress(upstream_meth, r_repeat_meth)
            x_line = np.array([min(upstream_meth), max(upstream_meth)])
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, 'r-', 
                   label=f'y={slope:.2f}x+{intercept:.2f} (r={r_value:.2f}, p={p_value:.4f})')
        
        ax.set_title("Upstream Coding Region vs r-repeat Methylation", fontsize=14)
        ax.set_xlabel("Upstream Coding Region Methylation (0-1)", fontsize=12)
        ax.set_ylabel("r-repeat Region Methylation (0-1)", fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(True, linestyle='--', alpha=0.7)
        if len(upstream_meth) > 1:  # Only add legend if regression line was added
            ax.legend()
        
        # Save figure
        upstream_meth_file = output_dir / "upstream_vs_r_repeat_methylation.png"
        plt.savefig(upstream_meth_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(upstream_meth_file))
        
        # Create another figure: Distribution of methylation differences
        fig, ax = plt.subplots(figsize=(10, 6))
        meth_diffs = [r - u for r, u in zip(r_repeat_meth, upstream_meth)]
        
        sns.histplot(meth_diffs, bins=30, kde=True, ax=ax)
        ax.axvline(x=0, color='r', linestyle='--', 
                  label='No difference between upstream and r-repeat')
        
        avg_diff = np.mean(meth_diffs)
        ax.axvline(x=avg_diff, color='g', linestyle='-', 
                  label=f'Mean difference: {avg_diff:.3f}')
        
        ax.set_title("Distribution of Methylation Differences (r-repeat - Upstream)", fontsize=14)
        ax.set_xlabel("Methylation Difference", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)
        ax.legend()
        
        # Save figure
        diff_file = output_dir / "methylation_difference_distribution.png"
        plt.savefig(diff_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(diff_file))
    
    # Add a new figure: Downstream IGS methylation vs r-repeat methylation
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Collect data for the plot
    downstream_meth = []
    r_repeat_meth = []
    for data in flattened_data:
        # Check if downstream IGS methylation data is available
        if ('downstream_igs_methylation' in data and 
            data['downstream_igs_methylation']['available']):
            downstream_avg = data['downstream_igs_methylation']['avg_methylation'] / 255  # Normalize to 0-1
            r_repeat_avg = data['avg_methylation'] / 255  # Normalize to 0-1
            
            # Only include points where downstream methylation could be calculated
            if downstream_avg >= 0:
                downstream_meth.append(downstream_avg)
                r_repeat_meth.append(r_repeat_avg)
    
    if downstream_meth and r_repeat_meth:
        # Scatter plot with regression line
        ax.scatter(downstream_meth, r_repeat_meth, alpha=0.6)
        
        # Add regression line
        if len(downstream_meth) > 1:  # Need at least 2 points for regression
            slope, intercept, r_value, p_value, std_err = linregress(downstream_meth, r_repeat_meth)
            x_line = np.array([min(downstream_meth), max(downstream_meth)])
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, 'r-', 
                   label=f'y={slope:.2f}x+{intercept:.2f} (r={r_value:.2f}, p={p_value:.4f})')
        
        ax.set_title("Downstream IGS Region vs r-repeat Methylation", fontsize=14)
        ax.set_xlabel("Downstream IGS Region Methylation (0-1)", fontsize=12)
        ax.set_ylabel("r-repeat Region Methylation (0-1)", fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(True, linestyle='--', alpha=0.7)
        if len(downstream_meth) > 1:  # Only add legend if regression line was added
            ax.legend()
        
        # Save figure
        downstream_meth_file = output_dir / "downstream_vs_r_repeat_methylation.png"
        plt.savefig(downstream_meth_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(downstream_meth_file))
        
    # Add a comparison figure: Upstream vs Downstream methylation
    if upstream_meth and downstream_meth and len(upstream_meth) == len(downstream_meth):
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Scatter plot with regression line
        ax.scatter(upstream_meth, downstream_meth, alpha=0.6)
        
        # Add regression line
        if len(upstream_meth) > 1:  # Need at least 2 points for regression
            slope, intercept, r_value, p_value, std_err = linregress(upstream_meth, downstream_meth)
            x_line = np.array([min(upstream_meth), max(upstream_meth)])
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, 'r-', 
                   label=f'y={slope:.2f}x+{intercept:.2f} (r={r_value:.2f}, p={p_value:.4f})')
        
        ax.set_title("Upstream Coding vs Downstream IGS Methylation", fontsize=14)
        ax.set_xlabel("Upstream Coding Region Methylation (0-1)", fontsize=12)
        ax.set_ylabel("Downstream IGS Region Methylation (0-1)", fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(True, linestyle='--', alpha=0.7)
        if len(upstream_meth) > 1:  # Only add legend if regression line was added
            ax.legend()
        
        # Save figure
        up_down_meth_file = output_dir / "upstream_vs_downstream_methylation.png"
        plt.savefig(up_down_meth_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(up_down_meth_file))
        
    # Add visualization that shows methylation levels across upstream, r-repeat, and downstream regions
    if upstream_meth and downstream_meth and r_repeat_meth:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Calculate average methylation levels
        avg_upstream = np.mean(upstream_meth)
        avg_r_repeat = np.mean(r_repeat_meth)
        avg_downstream = np.mean(downstream_meth)
        
        # Create bar plot
        regions = ['Upstream\nCoding', 'r-repeat', 'Downstream\nIGS']
        values = [avg_upstream, avg_r_repeat, avg_downstream]
        
        bars = ax.bar(regions, values, color=['skyblue', 'lightgreen', 'salmon'])
        
        # Add labels on top of bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                   f'{height:.3f}', ha='center', va='bottom')
        
        ax.set_title("Average Methylation Across Genomic Regions", fontsize=14)
        ax.set_xlabel("Region", fontsize=12)
        ax.set_ylabel("Average Methylation (0-1)", fontsize=12)
        ax.set_ylim(0, 1)
        ax.grid(axis='y', alpha=0.3)
        
        # Save figure
        regions_meth_file = output_dir / "methylation_across_regions.png"
        plt.savefig(regions_meth_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(regions_meth_file))
    
    # Add stats for downstream vs r-repeat methylation comparison
    if downstream_meth and r_repeat_meth:
        slope, intercept, r_value, p_value, std_err = linregress(downstream_meth, r_repeat_meth)
        
        visualization_stats["downstream_vs_r_repeat"] = {
            "correlation": {
                "r_value": float(r_value),
                "p_value": float(p_value)
            },
            "regression": {
                "slope": float(slope),
                "intercept": float(intercept),
                "std_err": float(std_err)
            },
            "mean_difference": float(np.mean([r - d for r, d in zip(r_repeat_meth, downstream_meth)])),
            "sample_count": len(downstream_meth)
        }
    
    # Save serialized statistics for cross-sample comparison
    with open(stats_dir / "methylation_visualization_stats.json", 'w') as f:
        json.dump(visualization_stats, f, indent=2)
    
    logger.info(f"Created {len(output_files)} methylation visualizations in {output_dir}")
    logger.info(f"Saved serialized visualization statistics to {stats_dir}")
    
    return output_files

def analyze_methylation(r_repeat_regions, mod_bam_path, output_dir=None, analyze_upstream=True, analyze_downstream=True):
    """
    Main function to analyze methylation in r-repeat regions, focusing on unit structure.
    
    Parameters:
        r_repeat_regions: Dictionary mapping read IDs to r-repeat regions
        mod_bam_path: Path to BAM file with modification (methylation) data
        output_dir: Directory for output files
        analyze_upstream: Whether to analyze upstream coding region methylation
        analyze_downstream: Whether to analyze downstream IGS region methylation
        
    Returns:
        tuple: (methylation data dict, list of visualization file paths, upstream analysis results, downstream analysis results)
    """
    # Set up proper directory structure
    mod_bam_path = Path(mod_bam_path)
    sample_name = mod_bam_path.stem
    
    if output_dir is None:
        output_dir = Path("r_repeat_results") / sample_name / "methylation"
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create data directory for serialized results
    data_dir = output_dir.parent / "data"
    data_dir.mkdir(exist_ok=True)
    
    # Extract methylation data
    r_repeat_data = extract_r_repeat_methylation(r_repeat_regions, mod_bam_path, output_dir)
    
    # Create visualizations
    viz_dir = output_dir / "figures"
    visualization_files = visualize_methylation_patterns(r_repeat_data, viz_dir)
    
    # Save primary methylation data to data directory
    with open(data_dir / "r_repeat_methylation.pkl", "wb") as f:
        pickle.dump(r_repeat_data, f)
    
    # Default results
    upstream_results = None
    downstream_results = None
    
    # Perform upstream coding methylation analysis if requested
    if analyze_upstream:
        upstream_dir = output_dir / "upstream_analysis"
        upstream_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("Performing upstream coding region methylation analysis...")
        upstream_results = analyze_upstream_coding_methylation(
            str(data_dir / "r_repeat_methylation.pkl"),
            output_dir=upstream_dir
        )
        
        if upstream_results and 'visualizations' in upstream_results:
            logger.info(f"Created {len(upstream_results['visualizations'])} upstream methylation visualizations")
    
    # Perform downstream IGS methylation analysis if requested
    if analyze_downstream:
        downstream_dir = output_dir / "downstream_analysis"
        downstream_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("Performing downstream IGS region methylation analysis...")
        downstream_results = analyze_downstream_igs_methylation(
            str(data_dir / "r_repeat_methylation.pkl"),
            output_dir=downstream_dir
        )
                
        if downstream_results and 'visualizations' in downstream_results:
            logger.info(f"Created {len(downstream_results['visualizations'])} downstream methylation visualizations")
    
    return r_repeat_data, visualization_files, upstream_results, downstream_results

if __name__ == "__main__":
    import argparse
    import json
    
    parser = argparse.ArgumentParser(description="Analyze methylation in r-repeat regions")
    parser.add_argument("--regions", "-r", required=True, help="Path to JSON file with r-repeat regions")
    parser.add_argument("--bam", "-b", required=True, help="Path to BAM file with methylation data")
    parser.add_argument("--output", "-o", help="Output directory for results")
    parser.add_argument("--upstream", "-u", action="store_true", help="Analyze upstream coding region methylation")
    parser.add_argument("--downstream", "-d", action="store_true", help="Analyze downstream IGS region methylation")
    
    args = parser.parse_args()
    
    # Load r-repeat regions
    with open(args.regions) as f:
        r_repeat_regions = json.load(f)
    
    analyze_methylation(r_repeat_regions, args.bam, args.output, args.upstream, args.downstream)
