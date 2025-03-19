"""
Multi-region analysis module for r-repeat regions.

This module analyzes patterns in reads that contain multiple r-repeat regions,
including relative positions, separation distances, and correlation between regions.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname=s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler("logs/r_repeat_multi_analysis.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def analyze_multi_regions(regions_file, output_dir=None):
    """
    Analyze reads with multiple r-repeat regions.
    
    Parameters:
        regions_file: Path to JSON file with r-repeat regions data
        output_dir: Directory for output files
        
    Returns:
        tuple: (multi-region data dictionary, visualization file paths dictionary)
    """
    # Set up output directory
    if output_dir is None:
        output_dir = Path("r_repeat_results") / "multi_region_analysis"
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load r-repeat regions data
    with open(regions_file, 'r') as f:
        all_regions = json.load(f)
    
    # Filter for reads with multiple r-repeat regions
    multi_regions = {read_id: regions for read_id, regions in all_regions.items()
                    if isinstance(regions, list) and len(regions) > 1}
    
    logger.info(f"Found {len(multi_regions)} reads with multiple r-repeat regions")
    
    if not multi_regions:
        logger.warning("No multi-region reads found. Analysis skipped.")
        return {}, {}
    
    # Analyze separation distances
    separations = []
    for read_id, regions in multi_regions.items():
        # Sort regions by start position
        sorted_regions = sorted(regions, key=lambda r: r['start'])
        
        # Calculate separations between adjacent regions
        for i in range(1, len(sorted_regions)):
            separation = sorted_regions[i]['start'] - sorted_regions[i-1]['end']
            if separation >= 0:  # Only consider non-overlapping regions
                separations.append(separation)
    
    # Create visualizations
    visualizations = {}
    
    # 1. Histogram of region separation distances
    if separations:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.histplot(separations, bins=30, kde=True, ax=ax)
        ax.set_title("Distribution of Separation Distances Between Adjacent r-repeat Regions", fontsize=14)
        ax.set_xlabel("Separation Distance (bp)", fontsize=12)
        ax.set_ylabel("Frequency", fontsize=12)
        
        # Add statistics
        stats_text = f"n = {len(separations)}\n"
        stats_text += f"min = {min(separations):.0f} bp\n"
        stats_text += f"mean = {np.mean(separations):.1f} bp\n"
        stats_text += f"median = {np.median(separations):.1f} bp\n"
        stats_text += f"max = {max(separations):.0f} bp"
        
        ax.annotate(stats_text, xy=(0.02, 0.95), xycoords='axes fraction',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
                    va='top', fontsize=10)
        
        # Save figure
        sep_file = output_dir / "separation_distances.png"
        plt.savefig(sep_file, dpi=300, bbox_inches="tight")
        plt.close()
        visualizations["separation_distances"] = str(sep_file)
    
    # 2. Bar plot of region count distribution
    region_counts = [len(regions) for regions in multi_regions.values()]
    count_dist = {}
    for count in region_counts:
        if count not in count_dist:
            count_dist[count] = 0
        count_dist[count] += 1
    
    if count_dist:
        fig, ax = plt.subplots(figsize=(10, 6))
        counts = sorted(count_dist.keys())
        freqs = [count_dist[count] for count in counts]
        
        bars = ax.bar(counts, freqs)
        
        # Add count labels on top of bars
        for bar, freq in zip(bars, freqs):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{freq}', ha='center', va='bottom')
        
        ax.set_title("Distribution of r-repeat Region Counts per Read", fontsize=14)
        ax.set_xlabel("Number of r-repeat Regions", fontsize=12)
        ax.set_ylabel("Number of Reads", fontsize=12)
        ax.set_xticks(counts)
        
        # Save figure
        count_file = output_dir / "region_count_distribution.png"
        plt.savefig(count_file, dpi=300, bbox_inches="tight")
        plt.close()
        visualizations["region_count_distribution"] = str(count_file)
    
    # Additional analysis and statistics
    stats = get_multi_region_statistics(multi_regions)
    
    # Save statistics to file
    stats_file = output_dir / "multi_region_statistics.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    logger.info(f"Multi-region analysis complete. Results saved to {output_dir}")
    return multi_regions, visualizations

def get_multi_region_statistics(multi_regions):
    """
    Calculate statistics for reads with multiple r-repeat regions.
    
    Parameters:
        multi_regions: Dictionary mapping read IDs to lists of r-repeat regions
        
    Returns:
        dict: Dictionary with statistics
    """
    if not multi_regions:
        return {"status": "no_data", "read_count": 0}
    
    # Count total regions
    region_count = sum(len(regions) for regions in multi_regions.values())
    
    # Count regions per read
    region_counts = [len(regions) for regions in multi_regions.values()]
    
    # Calculate separation distances
    separations = []
    overlaps = []
    for read_id, regions in multi_regions.items():
        # Sort regions by start position
        sorted_regions = sorted(regions, key=lambda r: r['start'])
        
        # Calculate separations between adjacent regions
        for i in range(1, len(sorted_regions)):
            separation = sorted_regions[i]['start'] - sorted_regions[i-1]['end']
            if separation >= 0:
                separations.append(separation)
            else:
                overlaps.append(-separation)  # Convert to positive overlap amount
    
    # Calculate region lengths
    region_lengths = []
    for regions in multi_regions.values():
        for region in regions:
            if 'length' in region:
                region_lengths.append(region['length'])
    
    # Calculate statistics for visualization and reporting
    stats = {
        "status": "success",
        "read_count": len(multi_regions),
        "region_count": region_count,
        "regions_per_read": {
            "min": min(region_counts),
            "max": max(region_counts),
            "mean": float(np.mean(region_counts)),
            "median": float(np.median(region_counts)),
            "distribution": {str(count): region_counts.count(count) for count in set(region_counts)}
        }
    }
    
    # Add separation statistics if available
    if separations:
        stats["separation_distances"] = {
            "count": len(separations),
            "min": float(min(separations)),
            "max": float(max(separations)),
            "mean": float(np.mean(separations)),
            "median": float(np.median(separations))
        }
    
    # Add overlap statistics if available
    if overlaps:
        stats["overlaps"] = {
            "count": len(overlaps),
            "min": float(min(overlaps)),
            "max": float(max(overlaps)),
            "mean": float(np.mean(overlaps)),
            "median": float(np.median(overlaps))
        }
    
    # Add region length statistics if available
    if region_lengths:
        stats["region_lengths"] = {
            "count": len(region_lengths),
            "min": float(min(region_lengths)),
            "max": float(max(region_lengths)),
            "mean": float(np.mean(region_lengths)),
            "median": float(np.median(region_lengths))
        }
    
    return stats

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze reads with multiple r-repeat regions")
    parser.add_argument("regions_file", help="Path to JSON file with r-repeat regions data")
    parser.add_argument("--output", "-o", help="Output directory for results")
    
    args = parser.parse_args()
    
    analyze_multi_regions(args.regions_file, args.output)
