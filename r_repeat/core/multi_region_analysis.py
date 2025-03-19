"""
Analysis of multiple r-repeat regions within single reads.

This module focuses on analyzing reads that contain multiple r-repeat regions,
investigating their relationship and spatial organization.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import json
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler("logs/r_repeat_multi_analysis.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def load_multi_region_data(r_repeat_regions_file):
    """
    Load r-repeat region data that contains multiple regions per read.
    
    Parameters:
        r_repeat_regions_file: Path to JSON file with r-repeat regions
        
    Returns:
        dict: Dictionary with multi-region r-repeat data
    """
    # Load r-repeat regions
    with open(r_repeat_regions_file) as f:
        r_repeat_regions = json.load(f)
    
    # Filter for reads with multiple regions
    multi_regions = {
        read_id: regions for read_id, regions in r_repeat_regions.items() 
        if isinstance(regions, list) and len(regions) > 1
    }
    
    logger.info(f"Found {len(multi_regions)}/{len(r_repeat_regions)} reads with multiple r-repeat regions")
    
    return multi_regions

def analyze_inter_repeat_distances(multi_regions):
    """
    Analyze distances between consecutive r-repeat regions in the same read.
    
    Parameters:
        multi_regions: Dictionary with multi-region r-repeat data
        
    Returns:
        list: List of distances between consecutive r-repeats
    """
    distances = []
    
    for read_id, regions in multi_regions.items():
        # Sort regions by start position
        sorted_regions = sorted(regions, key=lambda r: r['start'])
        
        # Calculate distances between consecutive regions
        for i in range(len(sorted_regions) - 1):
            current_end = sorted_regions[i]['end']
            next_start = sorted_regions[i+1]['start']
            distance = next_start - current_end
            distances.append(distance)
    
    logger.info(f"Analyzed {len(distances)} inter-repeat distances")
    logger.info(f"Distance statistics: min={min(distances) if distances else 'N/A'}, "
               f"mean={np.mean(distances) if distances else 'N/A':.1f}, "
               f"max={max(distances) if distances else 'N/A'}")
    
    return distances

def analyze_region_similarities(multi_regions):
    """
    Analyze similarities between r-repeat regions in the same read.
    
    Parameters:
        multi_regions: Dictionary with multi-region r-repeat data
        
    Returns:
        dict: Dictionary with similarity analysis results
    """
    length_diffs = []
    
    for read_id, regions in multi_regions.items():
        # Sort regions by start position
        sorted_regions = sorted(regions, key=lambda r: r['start'])
        
        # Calculate length differences between consecutive regions
        for i in range(len(sorted_regions) - 1):
            length1 = sorted_regions[i]['length']
            length2 = sorted_regions[i+1]['length']
            diff = abs(length1 - length2)
            length_diffs.append(diff)
    
    logger.info(f"Analyzed similarities between {len(length_diffs)} region pairs")
    
    return {
        'length_differences': length_diffs,
        'avg_length_diff': np.mean(length_diffs) if length_diffs else 0
    }

def get_multi_region_statistics(multi_regions):
    """
    Calculate statistics on multi-region r-repeats for serialization.
    
    Parameters:
        multi_regions: Dictionary with multi-region r-repeat data
        
    Returns:
        dict: Dictionary with statistics about multi-region patterns
    """
    if not multi_regions:
        return {}
    
    # Calculate distances between consecutive regions
    distances = []
    length_diffs = []
    
    for read_id, regions in multi_regions.items():
        # Sort regions by start position
        sorted_regions = sorted(regions, key=lambda r: r['start'])
        
        # Calculate distances and length differences
        for i in range(len(sorted_regions) - 1):
            current_end = sorted_regions[i]['end']
            next_start = sorted_regions[i+1]['start']
            distance = next_start - current_end
            distances.append(distance)
            
            length1 = sorted_regions[i]['length']
            length2 = sorted_regions[i+1]['length']
            diff = abs(length1 - length2)
            length_diffs.append(diff)
    
    # Calculate distribution of regions per read
    region_count_distribution = {}
    for read_id, regions in multi_regions.items():
        count = len(regions)
        if count not in region_count_distribution:
            region_count_distribution[count] = 0
        region_count_distribution[count] += 1
    
    # Convert to serializable format (keys must be strings in JSON)
    region_count_distribution = {
        str(count): freq 
        for count, freq in region_count_distribution.items()
    }
    
    # Compile statistics
    stats = {
        "region_count_distribution": region_count_distribution,
        "distances": {
            "count": len(distances),
            "min": min(distances) if distances else None,
            "mean": float(np.mean(distances)) if distances else None,
            "median": float(np.median(distances)) if distances else None,
            "max": max(distances) if distances else None,
            "std_dev": float(np.std(distances)) if distances else None
        },
        "length_differences": {
            "count": len(length_diffs),
            "min": min(length_diffs) if length_diffs else None,
            "mean": float(np.mean(length_diffs)) if length_diffs else None,
            "median": float(np.median(length_diffs)) if length_diffs else None,
            "max": max(length_diffs) if length_diffs else None,
            "std_dev": float(np.std(length_diffs)) if length_diffs else None
        }
    }
    
    return stats

def visualize_multi_region_analysis(multi_regions, distances, similarities, output_dir):
    """
    Create visualizations of multi-region r-repeat analysis.
    
    Parameters:
        multi_regions: Dictionary with multi-region r-repeat data
        distances: List of distances between consecutive r-repeats
        similarities: Dictionary with similarity analysis results
        output_dir: Directory for output files
        
    Returns:
        list: List of paths to generated visualization files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create stats directory for serialized visualization data
    stats_dir = output_dir.parent / "stats"
    stats_dir.mkdir(exist_ok=True)
    
    output_files = []
    visualization_stats = {}
    
    # Filter regions with problematic units for visualization purposes only
    filtered_multi_regions = {}
    for read_id, regions in multi_regions.items():
        # Only include the read if it contains at least one region without problematic units
        regions_without_problems = [r for r in regions if not r.get('has_problematic_unit', False)]
        if regions_without_problems:
            filtered_multi_regions[read_id] = regions_without_problems
    
    logger.info(f"Using {len(filtered_multi_regions)}/{len(multi_regions)} reads for visualization (excluding reads with only problematic units)")
    
    # Figure 1: Distribution of r-repeat counts per read
    fig, ax = plt.subplots(figsize=(10, 6))
    counts = [len(regions) for regions in filtered_multi_regions.values()]
    
    sns.histplot(counts, bins=range(2, max(counts)+2), kde=False, discrete=True, ax=ax)
    ax.set_title("Number of r-repeat Regions per Read", fontsize=14)
    ax.set_xlabel("r-repeat Count", fontsize=12)
    ax.set_ylabel("Number of Reads", fontsize=12)
    ax.set_xticks(range(2, max(counts)+1))
    
    counts_file = output_dir / "r_repeat_count_distribution.png"
    plt.savefig(counts_file, dpi=300, bbox_inches="tight")
    plt.close()
    output_files.append(str(counts_file))
    
    # Figure 2: Distribution of inter-repeat distances
    if distances:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.histplot(distances, bins=30, kde=True, ax=ax)
        ax.set_title("Distribution of Distances Between r-repeats", fontsize=14)
        ax.set_xlabel("Distance (bp)", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)
        
        distance_file = output_dir / "inter_repeat_distance_distribution.png"
        plt.savefig(distance_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(distance_file))
    
    # Figure 3: Distribution of length differences
    length_diffs = similarities['length_differences']
    if length_diffs:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.histplot(length_diffs, bins=30, kde=True, ax=ax)
        ax.set_title("Distribution of Length Differences Between Adjacent r-repeats", fontsize=14)
        ax.set_xlabel("Length Difference (bp)", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)
        
        diff_file = output_dir / "length_difference_distribution.png"
        plt.savefig(diff_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(diff_file))
    
    # Collect statistics for cross-sample comparison
    visualization_stats = {
        "region_counts": {
            "distribution": {str(count): counts.count(count) for count in set(counts)},
            "mean": float(np.mean(counts)) if counts else None,
            "median": float(np.median(counts)) if counts else None,
            "max": max(counts) if counts else None,
            "min": min(counts) if counts else None
        },
        "inter_repeat_distances": {
            "count": len(distances),
            "mean": float(np.mean(distances)) if distances else None,
            "median": float(np.median(distances)) if distances else None,
            "min": min(distances) if distances else None,
            "max": max(distances) if distances else None,
            "std_dev": float(np.std(distances)) if distances else None
        },
        "length_differences": {
            "count": len(similarities['length_differences']),
            "mean": float(np.mean(similarities['length_differences'])) if similarities['length_differences'] else None,
            "median": float(np.median(similarities['length_differences'])) if similarities['length_differences'] else None,
            "min": min(similarities['length_differences']) if similarities['length_differences'] else None,
            "max": max(similarities['length_differences']) if similarities['length_differences'] else None,
            "std_dev": float(np.std(similarities['length_differences'])) if similarities['length_differences'] else None
        },
        "sample_size": len(multi_regions)
    }
    
    # Save visualization statistics for cross-sample comparison
    with open(stats_dir / "multi_region_visualization_stats.json", 'w') as f:
        json.dump(visualization_stats, f, indent=2)
    
    logger.info(f"Created {len(output_files)} multi-region visualizations in {output_dir}")
    logger.info(f"Saved serialized visualization statistics to {stats_dir}")
    
    return output_files

def analyze_multi_regions(r_repeat_regions_file, output_dir=None):
    """
    Main function to analyze multiple r-repeat regions in reads.
    
    Parameters:
        r_repeat_regions_file: Path to JSON file with r-repeat regions
        output_dir: Directory for output files
        
    Returns:
        tuple: (multi-region data, list of visualization file paths)
    """
    # Set up proper directory structure
    regions_path = Path(r_repeat_regions_file)
    parent_dir = regions_path.parent
    
    # Try to extract sample name from parent directory structure
    sample_name = parent_dir.name
    if sample_name == "data":
        sample_name = parent_dir.parent.name
    
    if output_dir is None:
        output_dir = parent_dir.parent / "multi_region_analysis"
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create data directory for serialized results
    data_dir = parent_dir.parent / "data" if parent_dir.name == "data" else parent_dir
    data_dir.mkdir(exist_ok=True)
    
    # Load multi-region data
    multi_regions = load_multi_region_data(r_repeat_regions_file)
    
    # Skip analysis if no multi-region reads
    if not multi_regions:
        logger.warning("No reads with multiple r-repeat regions found")
        return multi_regions, []
    
    # Analyze inter-repeat distances
    distances = analyze_inter_repeat_distances(multi_regions)
    
    # Analyze region similarities
    similarities = analyze_region_similarities(multi_regions)
    
    # Create visualizations
    visualization_files = visualize_multi_region_analysis(
        multi_regions, distances, similarities, output_dir / "visualizations"
    )
    
    # Save comprehensive results to data directory
    results = {
        "sample_name": sample_name,
        "multi_regions": multi_regions,
        "statistics": get_multi_region_statistics(multi_regions)
    }
    
    output_file = data_dir / "multi_region_analysis.json"
    with open(output_file, "w") as f:
        # Convert multi_regions data to serializable format
        serializable_results = {
            "sample_name": sample_name,
            "multi_region_count": len(multi_regions),
            "statistics": results["statistics"]
        }
        json.dump(serializable_results, f, indent=2)
    
    # Save full data with pickle
    with open(data_dir / "multi_region_full_data.pkl", "wb") as f:
        pickle.dump(results, f)
    
    logger.info(f"Multi-region analysis complete. Results saved to {output_file}")
    
    return multi_regions, visualization_files

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze multiple r-repeat regions in reads")
    parser.add_argument("--regions", "-r", required=True, help="Path to JSON file with r-repeat regions")
    parser.add_argument("--output", "-o", help="Output directory for results")
    
    args = parser.parse_args()
    
    analyze_multi_regions(args.regions, args.output)
