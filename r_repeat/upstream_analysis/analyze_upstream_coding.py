"""
Analysis module for upstream coding region methylation patterns in relation to r-repeat regions.

This module is designed to be used as part of the r-repeat methylation analysis pipeline
but can also be run independently for more detailed upstream coding region analysis.
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import argparse
import logging
import json  # Add the missing json import

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

def analyze_upstream_coding_methylation(data_file, output_dir=None, low_upstream_threshold=0.3):
    """
    Analyze the relationship between upstream coding region methylation and r-repeat methylation.
    
    Parameters:
        data_file: Path to a pickle file containing r-repeat methylation data with upstream information
        output_dir: Directory for output files
        low_upstream_threshold: Threshold to define low upstream methylation (default: 0.3)
        
    Returns:
        dict: Summary statistics of the analysis
    """
    # Set up proper directory structure
    data_path = Path(data_file)
    parent_dir = data_path.parent
    
    # Try to extract sample name from file path
    sample_name = parent_dir.name
    if sample_name == "data":
        sample_name = parent_dir.parent.name
    
    if output_dir is None:
        output_dir = parent_dir.parent / "upstream_analysis"
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create stats directory for serialized comparison data
    stats_dir = output_dir / "stats"
    stats_dir.mkdir(exist_ok=True)
    
    # Load the data
    logger.info(f"Loading data from {data_file}")
    with open(data_file, 'rb') as f:
        r_repeat_data = pickle.load(f)
    
    # Check if we have multiple regions per read
    has_multiple_regions = any(isinstance(data, list) for data in r_repeat_data.values())
    
    # Extract methylation data pairs
    upstream_meth = []
    r_repeat_meth = []
    region_count = 0
    
    # New lists to store regions with low upstream methylation
    low_upstream_regions = []
    
    if has_multiple_regions:
        # Handle multiple regions per read
        for read_id, regions in r_repeat_data.items():
            for region in regions:
                if ('upstream_coding_methylation' in region and 
                    region['upstream_coding_methylation']['available']):
                    
                    upstream_avg = region['upstream_coding_methylation']['avg_methylation'] / 255
                    r_repeat_avg = region['avg_methylation'] / 255
                    
                    if upstream_avg >= 0:
                        upstream_meth.append(upstream_avg)
                        r_repeat_meth.append(r_repeat_avg)
                        region_count += 1
                        
                        # Store regions with low upstream methylation
                        if upstream_avg < low_upstream_threshold:
                            low_upstream_regions.append({
                                'read_id': read_id,
                                'region_data': region,
                                'upstream_methylation': upstream_avg,
                                'r_repeat_methylation': r_repeat_avg
                            })
    else:
        # Handle single region per read
        for read_id, data in r_repeat_data.items():
            if ('upstream_coding_methylation' in data and 
                data['upstream_coding_methylation']['available']):
                
                upstream_avg = data['upstream_coding_methylation']['avg_methylation'] / 255
                r_repeat_avg = data['avg_methylation'] / 255
                
                if upstream_avg >= 0:
                    upstream_meth.append(upstream_avg)
                    r_repeat_meth.append(r_repeat_avg)
                    region_count += 1
                    
                    # Store regions with low upstream methylation
                    if upstream_avg < low_upstream_threshold:
                        low_upstream_regions.append({
                            'read_id': read_id,
                            'region_data': data,
                            'upstream_methylation': upstream_avg,
                            'r_repeat_methylation': r_repeat_avg
                        })
    
    logger.info(f"Analyzed {region_count} regions with available upstream coding data")
    logger.info(f"Found {len(low_upstream_regions)} regions with low upstream methylation (<{low_upstream_threshold})")
    
    if len(upstream_meth) < 2:
        logger.warning("Not enough data points with upstream coding information for analysis")
        return {
            'status': 'error',
            'message': 'Not enough data points with upstream coding information'
        }
    
    # Calculate basic statistics
    upstream_avg = np.mean(upstream_meth)
    r_repeat_avg = np.mean(r_repeat_meth)
    meth_diffs = [r - u for r, u in zip(r_repeat_meth, upstream_meth)]
    avg_diff = np.mean(meth_diffs)
    
    # Calculate correlation
    from scipy.stats import pearsonr, spearmanr
    pearson_r, pearson_p = pearsonr(upstream_meth, r_repeat_meth)
    spearman_r, spearman_p = spearmanr(upstream_meth, r_repeat_meth)
    
    # Perform linear regression
    from scipy.stats import linregress
    slope, intercept, r_value, p_value, std_err = linregress(upstream_meth, r_repeat_meth)
    
    # Plot 1: Scatter plot with regression line
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(upstream_meth, r_repeat_meth, alpha=0.6)
    
    # Add regression line
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
    ax.legend()
    
    scatter_file = output_dir / "upstream_vs_r_repeat_scatter.png"
    plt.savefig(scatter_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Plot 2: Distribution of methylation differences
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(meth_diffs, bins=30, kde=True, ax=ax)
    
    ax.axvline(x=0, color='r', linestyle='--', 
              label='No difference between upstream and r-repeat')
    ax.axvline(x=avg_diff, color='g', linestyle='-', 
              label=f'Mean difference: {avg_diff:.3f}')
    
    ax.set_title("Distribution of Methylation Differences (r-repeat - Upstream)", fontsize=14)
    ax.set_xlabel("Methylation Difference", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.legend()
    
    diff_file = output_dir / "methylation_difference_distribution.png"
    plt.savefig(diff_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Plot 3: Heatmap of methylation status (binned)
    # Convert continuous methylation levels to binned categories for visualization
    df = pd.DataFrame({
        'Upstream': upstream_meth,
        'r-repeat': r_repeat_meth
    })
    
    # Create bins for methylation levels
    bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    labels = ['0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0']
    
    df['Upstream_binned'] = pd.cut(df['Upstream'], bins=bins, labels=labels)
    df['r-repeat_binned'] = pd.cut(df['r-repeat'], bins=bins, labels=labels)
    
    # Create contingency table
    ctab = pd.crosstab(df['Upstream_binned'], df['r-repeat_binned'], normalize='all') * 100
    
    # Plot heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(ctab, annot=True, cmap="YlGnBu", fmt=".1f", cbar_kws={'label': 'Percentage (%)'})
    
    ax.set_title("Relationship Between Upstream and r-repeat Methylation Levels", fontsize=14)
    ax.set_xlabel("r-repeat Region Methylation", fontsize=12)
    ax.set_ylabel("Upstream Coding Region Methylation", fontsize=12)
    
    heatmap_file = output_dir / "methylation_relationship_heatmap.png"
    plt.savefig(heatmap_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Create summary statistics
    summary = {
        'status': 'success',
        'upstream_avg_methylation': upstream_avg,
        'r_repeat_avg_methylation': r_repeat_avg,
        'mean_difference': avg_diff,
        'pearson_correlation': {
            'r': pearson_r,
            'p_value': pearson_p
        },
        'spearman_correlation': {
            'rho': spearman_r,
            'p_value': spearman_p
        },
        'linear_regression': {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_value**2,
            'p_value': p_value,
            'std_error': std_err
        },
        'region_count': region_count,
        'visualizations': {
            'scatter_plot': str(scatter_file),
            'difference_histogram': str(diff_file),
            'relationship_heatmap': str(heatmap_file)
        }
    }
    
    # New: Create visualizations for regions with low upstream methylation
    low_upstream_visualizations = visualize_low_upstream_regions(
        low_upstream_regions, 
        output_dir, 
        low_upstream_threshold
    )
    
    # Add new visualizations to summary
    summary['visualizations'].update(low_upstream_visualizations)
    
    # New: Create visualization for upstream methylation distribution by unit count
    unit_count_visualization = visualize_upstream_methylation_by_unit_count(
        r_repeat_data, 
        output_dir, 
        has_multiple_regions
    )
    
    if unit_count_visualization:
        summary['visualizations']['upstream_by_unit_count'] = unit_count_visualization
    
    # Create comprehensive serialized data for cross-sample comparison
    comparison_data = {
        "sample_name": sample_name,
        "region_count": region_count,
        "low_upstream_count": len(low_upstream_regions),
        "low_upstream_threshold": low_upstream_threshold,
        "correlations": {
            "pearson": {
                "r": float(pearson_r),
                "p_value": float(pearson_p)
            },
            "spearman": {
                "rho": float(spearman_r),
                "p_value": float(spearman_p)
            }
        },
        "regression": {
            "slope": float(slope),
            "intercept": float(intercept),
            "r_squared": float(r_value**2),
            "p_value": float(p_value)
        },
        "methylation_levels": {
            "upstream": {
                "mean": float(upstream_avg),
                "quartiles": {
                    "q1": float(np.percentile(upstream_meth, 25)),
                    "median": float(np.median(upstream_meth)),
                    "q3": float(np.percentile(upstream_meth, 75))
                }
            },
            "r_repeat": {
                "mean": float(r_repeat_avg),
                "quartiles": {
                    "q1": float(np.percentile(r_repeat_meth, 25)),
                    "median": float(np.median(r_repeat_meth)),
                    "q3": float(np.percentile(r_repeat_meth, 75))
                }
            },
            "mean_difference": float(avg_diff)
        },
        "unit_count_distribution": {}
    }
    
    # Add unit count distribution if available
    unit_count_data = {}
    for region_info in low_upstream_regions:
        region_data = region_info['region_data']
        if 'coordinates' in region_data and 'unit_count' in region_data['coordinates']:
            unit_count = region_data['coordinates']['unit_count']
            if unit_count not in unit_count_data:
                unit_count_data[unit_count] = 0
            unit_count_data[unit_count] += 1
    
    comparison_data["unit_count_distribution"] = {str(k): v for k, v in unit_count_data.items()}
    
    # Save comprehensive comparison data
    comparison_file = stats_dir / "upstream_comparison_data.json"
    with open(comparison_file, 'w') as f:
        json.dump(comparison_data, f, indent=2)
    
    # Save summary to file
    summary_file = output_dir / "upstream_analysis_summary.json"
    with open(summary_file, 'w') as f:
        # Convert numpy values to standard Python types for JSON serialization
        json_safe_summary = {
            k: ({kk: float(vv) if isinstance(vv, (np.float64, np.float32)) else vv 
                for kk, vv in v.items()} if isinstance(v, dict) else 
                (float(v) if isinstance(v, (np.float64, np.float32)) else v))
            for k, v in summary.items()
        }
        json.dump(json_safe_summary, f, indent=2)
    
    logger.info(f"Analysis complete. Results saved to {output_dir}")
    logger.info(f"Correlation: Pearson r={pearson_r:.3f} (p={pearson_p:.4f}), " 
               f"Spearman rho={spearman_r:.3f} (p={spearman_p:.4f})")
    logger.info(f"Mean difference (r-repeat - upstream): {avg_diff:.3f}")
    logger.info(f"Cross-sample comparison data saved to {comparison_file}")
    
    # Add unit-specific methylation analysis
    unit_analysis_results = analyze_unit_specific_methylation(
        r_repeat_data, 
        output_dir, 
        has_multiple_regions
    )
    
    if unit_analysis_results and unit_analysis_results['status'] == 'success':
        summary['unit_methylation'] = unit_analysis_results['unit_stats']
        # Add visualizations to summary
        summary['visualizations'].update(unit_analysis_results['visualizations'])
    
    return summary

def visualize_low_upstream_regions(low_upstream_regions, output_dir, threshold):
    """
    Create visualizations specifically for regions with low upstream methylation.
    
    Parameters:
        low_upstream_regions: List of regions with low upstream methylation
        output_dir: Directory for output files
        threshold: Threshold used to define low upstream methylation
        
    Returns:
        dict: Dictionary of paths to visualization files
    """
    visualizations = {}
    
    if not low_upstream_regions:
        logger.warning("No regions with low upstream methylation to visualize")
        return visualizations
    
    # Organize regions by unit count
    regions_by_unit_count = {}
    for region_info in low_upstream_regions:
        region_data = region_info['region_data']
        
        # Skip regions with no unit count information
        if 'coordinates' not in region_data or 'unit_count' not in region_data['coordinates']:
            continue
            
        unit_count = region_data['coordinates']['unit_count']
        
        if unit_count not in regions_by_unit_count:
            regions_by_unit_count[unit_count] = []
        
        regions_by_unit_count[unit_count].append(region_data)
    
    # Visualize methylation patterns along r-repeat by unit count
    # (only for regions with low upstream methylation)
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for unit_count, regions in sorted(regions_by_unit_count.items()):
        if unit_count < 1 or len(regions) < 3:  # Skip if no meaningful unit count or too few samples
            continue
            
        # Group methylation data by relative position
        all_positions = {}
        for region in regions:
            if 'methylation_data' not in region:
                continue
                
            length = region['coordinates']['length']
            for pos, score in region['methylation_data']:
                rel_pos = int((pos / length) * 100)  # Convert to 0-100% position
                if rel_pos not in all_positions:
                    all_positions[rel_pos] = []
                all_positions[rel_pos].append(score / 255)  # Normalize score to 0-1
        
        if not all_positions:
            continue
            
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
    
    ax.set_title(f"Methylation Patterns in r-repeat Regions with Low Upstream Methylation (<{threshold})", 
                fontsize=14)
    ax.set_xlabel("Relative Position in r-repeat Region (%)", fontsize=12)
    ax.set_ylabel("Average Methylation (0-1)", fontsize=12)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(title="Unit count", loc="best")
    
    low_pattern_file = output_dir / "low_upstream_methylation_pattern.png"
    plt.savefig(low_pattern_file, dpi=300, bbox_inches="tight")
    plt.close()
    visualizations['low_upstream_pattern'] = str(low_pattern_file)
    
    # NEW: Create unit-based methylation visualization
    fig, axs = plt.subplots(figsize=(16, 12), nrows=min(5, len(regions_by_unit_count)), ncols=1, 
                           constrained_layout=True, sharex=True)
    
    # If there's only one subplot, convert to list for consistent indexing
    if len(regions_by_unit_count) == 1:
        axs = [axs]
    
    # Sort unit counts
    sorted_unit_counts = sorted(regions_by_unit_count.keys())
    
    # Limit to at most 5 unit counts for readability
    unit_counts_to_plot = sorted_unit_counts[:5]
    
    # Create color palette for the plots
    colors = plt.cm.viridis(np.linspace(0, 1, max(unit_counts_to_plot)+1))
    
    for i, unit_count in enumerate(unit_counts_to_plot):
        if unit_count < 1 or len(regions_by_unit_count[unit_count]) < 3:
            continue
            
        ax = axs[i]
        regions = regions_by_unit_count[unit_count]
        
        # Extract methylation data for each unit separately
        unit_methylation = {k: [] for k in range(1, unit_count+1)}
        unit_methylation['upstream'] = []
        
        # First collect raw methylation data from each region
        for region in regions:
            # Process r-repeat methylation data (if unit-specific data not available)
            # This ensures we at least have data to show even when unit_methylation is missing
            if 'methylation_data' in region and 'coordinates' in region and 'length' in region['coordinates']:
                length = region['coordinates']['length']
                # Approximate unit length by dividing by unit count
                approx_unit_length = length / unit_count if unit_count > 0 else length
                
                for pos, score in region['methylation_data']:
                    # Determine which unit this position belongs to
                    unit_num = min(int(pos / approx_unit_length) + 1, unit_count)
                    
                    # Calculate position relative to the start of this unit
                    rel_unit_pos = pos - (unit_num - 1) * approx_unit_length
                    
                    # Normalize to 0-100% within the unit
                    rel_pos = int((rel_unit_pos / approx_unit_length) * 100)
                    
                    # Store the methylation data in the appropriate unit
                    if unit_num in unit_methylation:
                        unit_methylation[unit_num].append({
                            'position': rel_pos,
                            'methylation': score / 255  # Normalize to 0-1
                        })
            
            # Extract unit-specific methylation (if available)
            if 'unit_methylation' in region and region['unit_methylation']:
                for unit_key, unit_data in region['unit_methylation'].items():
                    # Extract unit number from key (e.g., "unit_1" -> 1)
                    unit_num = int(unit_key.split('_')[1])
                    
                    if unit_num <= unit_count:
                        # Normalize positions to 0-100% within the unit
                        unit_positions = {}
                        unit_length = unit_data.get('unit_length', 0)
                        
                        if unit_length > 0:
                            for pos, score in unit_data.get('methylation_data', []):
                                rel_pos = int((pos / unit_length) * 100)
                                if rel_pos not in unit_positions:
                                    unit_positions[rel_pos] = []
                                unit_positions[rel_pos].append(score / 255)  # Normalize to 0-1
                                
                            # Calculate average methylation at each position
                            for pos, scores in unit_positions.items():
                                unit_methylation[unit_num].append({
                                    'position': pos,
                                    'methylation': np.mean(scores)
                                })
            
            # Add upstream methylation data
            if ('upstream_coding_methylation' in region and 
                region['upstream_coding_methylation']['available'] and
                region['upstream_coding_methylation']['methylation_data']):
                
                upstream_positions = {}
                upstream_length = 2000  # Standard upstream analysis length
                
                for pos, score in region['upstream_coding_methylation']['methylation_data']:
                    rel_pos = int((pos / upstream_length) * 100)
                    if rel_pos not in upstream_positions:
                        upstream_positions[rel_pos] = []
                    upstream_positions[rel_pos].append(score / 255)
                
                # Calculate average methylation at each position
                for pos, scores in upstream_positions.items():
                    unit_methylation['upstream'].append({
                        'position': pos,
                        'methylation': np.mean(scores)
                    })
        
        # Add shaded rectangles to indicate the motif regions in each unit
        # Use different motif lengths based on unit position and count
        for unit_num in range(1, unit_count+1):
            # Determine motif length ratio based on unit and total units
            if unit_count == 3 and unit_num == 3:
                # Third unit in 3-unit regions is typically shorter
                motif_length_ratio = 350 / 650  # Shorter motif
            elif unit_num == unit_count and unit_count > 1:
                # Last unit is often slightly shorter
                motif_length_ratio = 380 / 650
            else:
                # Standard motif length for most units
                motif_length_ratio = 400 / 650
            
            # Calculate the motif region coordinates (0-100 scale per unit)
            motif_start = (unit_num - 1) * 100
            motif_end = motif_start + motif_length_ratio * 100
            
            # Add shaded rectangle for the motif
            ax.axvspan(motif_start, motif_end, alpha=0.15, color='blue', 
                      label='Motif region' if unit_num == 1 else "")
            
            # Add CT-region annotation
            ct_region_start = motif_end
            ct_region_end = unit_num * 100
            ax.axvspan(ct_region_start, ct_region_end, alpha=0.15, color='green',
                      label='CT-repeat region' if unit_num == 1 else "")
        
        # Plot methylation pattern for each unit
        for unit_num in range(1, unit_count+1):
            if unit_num in unit_methylation and unit_methylation[unit_num]:
                # Convert to DataFrame for easier plotting
                df = pd.DataFrame(unit_methylation[unit_num])
                
                # Check if we have enough data points to bin
                if len(df) > 0:
                    # Bin positions to reduce noise - adjust number of bins based on data
                    num_bins = min(20, max(5, len(df) // 3))  # At least 5 bins, at most 20
                    df['position_bin'] = pd.cut(df['position'], bins=num_bins, labels=False)
                    
                    # Group by position bin and calculate mean
                    binned_data = df.groupby('position_bin', observed=True)['methylation'].mean().reset_index()
                    
                    # Plot only if we have data after binning
                    if not binned_data.empty:
                        # Calculate x-positions based on unit position
                        # Shift each unit to its proper position on the x-axis
                        x_positions = (unit_num - 1) * 100 + binned_data['position_bin'] * (100.0 / num_bins)
                        
                        ax.plot(x_positions, binned_data['methylation'], 
                               label=f"Unit {unit_num}", 
                               color=colors[unit_num], 
                               linewidth=2, alpha=0.8)
        
        # Add upstream data if available
        if unit_methylation['upstream']:
            df = pd.DataFrame(unit_methylation['upstream'])
            df['position_bin'] = pd.cut(df['position'], bins=20, labels=False)
            binned_data = df.groupby('position_bin')['methylation'].mean().reset_index()
            
            # Fix the upstream data plotting - create a proper x-axis array
            # Convert from upstream position bins to plot coordinates
            x_positions = -100 + binned_data['position_bin'] * 5  # Convert to negative positions
            
            # Plot upstream data with dashed line - use the correct array shapes
            ax.plot(x_positions, binned_data['methylation'], 
                   label="Upstream", 
                   color='gray', 
                   linestyle='--',
                   linewidth=2, alpha=0.8)
        
        # Add vertical line separating upstream and r-repeat
        ax.axvline(x=0, color='black', linestyle='-', alpha=0.5, label='r-repeat start')
        
        # Add annotations for unit boundaries based on typical lengths
        for u in range(1, unit_count):
            ax.axvline(x=u*100, color='gray', linestyle=':', alpha=0.3)
        
        ax.set_title(f"{unit_count} units (n={len(regions)})", fontsize=12)
        ax.set_ylim(0, 1)
        ax.set_xlim(-100, 100*unit_count)  # Set proper x-axis limits to show all units
        ax.grid(True, linestyle='--', alpha=0.3)
        
        # Add informative text labels for units
        for u in range(1, unit_count+1):
            ax.text((u-0.5)*100, 0.95, f"Unit {u}", ha='center', fontsize=9,
                   bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray', boxstyle='round,pad=0.2'))
        
        # Only add x-axis label to bottom subplot
        if i == len(axs)-1:
            ax.set_xlabel("Relative Position (%)", fontsize=12)
        
        # Add y-axis label to all subplots
        ax.set_ylabel("Average Methylation (0-1)", fontsize=10)
        
        # Add legend to each subplot, positioned inside the plot area without overlaying axis labels
        # Position in the upper left corner of the actual plot area
        ax.legend(title="Region", loc="upper left", fontsize=9, 
                 bbox_to_anchor=(0.01, 0.99), framealpha=0.9)
    
    # Adjust the figure layout to make room for the legends
    fig.tight_layout()
    
    # Add overall title
    fig.suptitle(f"Methylation Patterns by Unit in r-repeats with Low Upstream Methylation (<{threshold})", fontsize=14)
    
    # Additional padding to accommodate the title
    plt.subplots_adjust(top=0.95)
    
    # Save figure
    unit_pattern_file = output_dir / "low_upstream_methylation_by_unit.png"
    plt.savefig(unit_pattern_file, dpi=300, bbox_inches="tight")
    plt.close()
    visualizations['low_upstream_by_unit'] = str(unit_pattern_file)
    
    # Create statistical summary for low upstream methylation regions
    if low_upstream_regions:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Extract r-repeat methylation values from low upstream regions
        low_upstream_r_repeat_meth = [r['r_repeat_methylation'] for r in low_upstream_regions]
        
        # Create histogram
        sns.histplot(low_upstream_r_repeat_meth, bins=20, kde=True, ax=ax)
        ax.axvline(x=np.mean(low_upstream_r_repeat_meth), color='r', linestyle='--',
                  label=f'Mean: {np.mean(low_upstream_r_repeat_meth):.3f}')
        
        ax.set_title(f"R-repeat Methylation Distribution (Upstream Methylation <{threshold})", 
                    fontsize=14)
        ax.set_xlabel("R-repeat Methylation Level (0-1)", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)
        ax.legend()
        
        meth_dist_file = output_dir / "low_upstream_r_repeat_methylation_distribution.png"
        plt.savefig(meth_dist_file, dpi=300, bbox_inches="tight")
        plt.close()
        visualizations['low_upstream_meth_distribution'] = str(meth_dist_file)
    
    logger.info(f"Created {len(visualizations)} visualizations for regions with low upstream methylation")
    return visualizations

def visualize_upstream_methylation_by_unit_count(r_repeat_data, output_dir, has_multiple_regions):
    """
    Create visualizations showing upstream methylation distribution by unit count.
    
    Parameters:
        r_repeat_data: Dictionary with r-repeat methylation data
        output_dir: Directory for output files
        has_multiple_regions: Whether the data has multiple regions per read
        
    Returns:
        str: Path to visualization file, or None if visualization cannot be created
    """
    # Collect upstream methylation data by unit count
    unit_count_upstream_meth = {}
    
    if has_multiple_regions:
        # Handle multiple regions per read
        for read_id, regions in r_repeat_data.items():
            for region in regions:
                if ('upstream_coding_methylation' in region and 
                    region['upstream_coding_methylation']['available'] and
                    'coordinates' in region and 
                    'unit_count' in region['coordinates']):
                    
                    unit_count = region['coordinates']['unit_count']
                    if unit_count < 1:  # Skip regions with no units
                        continue
                        
                    upstream_avg = region['upstream_coding_methylation']['avg_methylation'] / 255
                    
                    if upstream_avg >= 0:  # Skip unavailable upstream data
                        if unit_count not in unit_count_upstream_meth:
                            unit_count_upstream_meth[unit_count] = []
                        unit_count_upstream_meth[unit_count].append(upstream_avg)
    else:
        # Handle single region per read
        for read_id, data in r_repeat_data.items():
            if ('upstream_coding_methylation' in data and 
                data['upstream_coding_methylation']['available'] and
                'coordinates' in data and 
                'unit_count' in data['coordinates']):
                
                unit_count = data['coordinates']['unit_count']
                if unit_count < 1:  # Skip regions with no units
                    continue
                    
                upstream_avg = data['upstream_coding_methylation']['avg_methylation'] / 255
                
                if upstream_avg >= 0:  # Skip unavailable upstream data
                    if unit_count not in unit_count_upstream_meth:
                        unit_count_upstream_meth[unit_count] = []
                    unit_count_upstream_meth[unit_count].append(upstream_avg)
    
    # Check if we have enough data
    if not unit_count_upstream_meth:
        logger.warning("No upstream methylation data by unit count available")
        return None
    
    # Remove unit counts with too few samples for meaningful visualization
    unit_count_upstream_meth = {k: v for k, v in unit_count_upstream_meth.items() if len(v) >= 3}
    
    if not unit_count_upstream_meth:
        logger.warning("Not enough samples per unit count for meaningful visualization")
        return None
    
    # Create boxplot visualization
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Prepare data for boxplot
    data = []
    for unit_count, values in sorted(unit_count_upstream_meth.items()):
        for val in values:
            data.append({
                'Unit Count': unit_count,
                'Upstream Methylation': val
            })
    
    df = pd.DataFrame(data)
    
    # Create boxplot
    sns.boxplot(x='Unit Count', y='Upstream Methylation', data=df, ax=ax)
    
    # Add sample sizes to x-axis labels
    labels = [f"{count}\n(n={len(unit_count_upstream_meth[count])})" 
              for count in sorted(unit_count_upstream_meth.keys())]
    ax.set_xticklabels(labels)
    
    ax.set_title("Upstream Coding Region Methylation by r-repeat Unit Count", fontsize=14)
    ax.set_xlabel("Number of r-repeat Units", fontsize=12)
    ax.set_ylabel("Upstream Region Methylation (0-1)", fontsize=12)
    ax.set_ylim(0, 1)
    
    # Add horizontal line at low methylation threshold
    ax.axhline(y=0.3, color='r', linestyle='--', label='Low methylation threshold (0.3)')
    ax.legend()
    
    # Save figure
    output_file = output_dir / "upstream_methylation_by_unit_count.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"Created upstream methylation by unit count visualization: {output_file}")
    return str(output_file)

def analyze_unit_specific_methylation(r_repeat_data, output_dir, has_multiple_regions):
    """
    Analyze methylation patterns specifically at the unit level for r-repeat regions.
    
    Parameters:
        r_repeat_data: Dictionary with r-repeat methylation data
        output_dir: Directory for output files
        has_multiple_regions: Whether the data has multiple regions per read
        
    Returns:
        dict: Dictionary with unit methylation analysis results and visualization paths
    """
    # Create stats directory for serialized comparison data
    stats_dir = output_dir / "stats"
    stats_dir.mkdir(exist_ok=True)
    
    # Dictionary to track unit-level methylation by unit count
    unit_methylation_by_count = {}
    
    # Extract unit methylation data from r-repeat data
    if has_multiple_regions:
        # Handle multiple regions per read
        for read_id, regions in r_repeat_data.items():
            for region in regions:
                if 'unit_methylation' in region and region['unit_methylation']:  # Check if unit_methylation exists and is not None
                    unit_count = region['coordinates'].get('unit_count', 0)
                    
                    if unit_count < 1:
                        continue
                        
                    if unit_count not in unit_methylation_by_count:
                        unit_methylation_by_count[unit_count] = {i+1: [] for i in range(unit_count)}
                    
                    for unit_id, unit_data in region['unit_methylation'].items():
                        unit_num = int(unit_id.split('_')[1])
                        if unit_num <= unit_count:
                            avg_meth = unit_data.get('avg_methylation', 0) / 255  # Normalize to 0-1
                            unit_methylation_by_count[unit_count][unit_num].append(avg_meth)
    else:
        # Handle single region per read
        for read_id, data in r_repeat_data.items():
            if 'unit_methylation' in data and data['unit_methylation']:  # Check if unit_methylation exists and is not None
                unit_count = data['coordinates'].get('unit_count', 0)
                
                if unit_count < 1:
                    continue
                    
                if unit_count not in unit_methylation_by_count:
                    unit_methylation_by_count[unit_count] = {i+1: [] for i in range(unit_count)}
                
                for unit_id, unit_data in data['unit_methylation'].items():
                    unit_num = int(unit_id.split('_')[1])
                    if unit_num <= unit_count:
                        avg_meth = unit_data.get('avg_methylation', 0) / 255  # Normalize to 0-1
                        unit_methylation_by_count[unit_count][unit_num].append(avg_meth)
    
    # Filter to include only unit counts with sufficient data
    unit_methylation_by_count = {
        count: unit_data for count, unit_data in unit_methylation_by_count.items()
        if all(len(values) >= 3 for values in unit_data.values())  # Require at least 3 samples per unit
    }
    
    if not unit_methylation_by_count:
        logger.warning("Not enough unit methylation data for analysis")
        return {'status': 'warning', 'message': 'Not enough unit methylation data'}
    
    # Create a figure to visualize methylation patterns by unit position for different unit counts
    visualizations = {}
    
    # Sort unit counts and limit to a reasonable number for visualization
    sorted_counts = sorted(unit_methylation_by_count.keys())
    counts_to_plot = sorted_counts[:min(6, len(sorted_counts))]
    
    fig, axs = plt.subplots(len(counts_to_plot), 1, figsize=(12, 4*len(counts_to_plot)), 
                           constrained_layout=True, sharex=True)
    
    # Handle single subplot case
    if len(counts_to_plot) == 1:
        axs = [axs]
    
    # Create boxplot for each unit count
    for i, unit_count in enumerate(counts_to_plot):
        ax = axs[i]
        unit_data = unit_methylation_by_count[unit_count]
        
        # Prepare data for boxplot
        plot_data = []
        for unit_num, values in unit_data.items():
            for val in values:
                plot_data.append({
                    'Unit Position': f"Unit {unit_num}",
                    'Methylation': val
                })
        
        df = pd.DataFrame(plot_data)
        
        # Create boxplot
        sns.boxplot(x='Unit Position', y='Methylation', data=df, ax=ax)
        
        # Add sample size to title
        sample_sizes = [len(unit_data[unit]) for unit in range(1, unit_count+1)]
        ax.set_title(f"{unit_count} units (n={min(sample_sizes)} to {max(sample_sizes)} regions per unit)")
        
        # Only add x-label to bottom subplot
        if i == len(axs)-1:
            ax.set_xlabel("Unit Position", fontsize=12)
        else:
            ax.set_xlabel("")
            
        ax.set_ylabel("Methylation (0-1)", fontsize=10)
        ax.set_ylim(0, 1)
    
    # Add overall title
    fig.suptitle("r-repeat Methylation by Unit Position", fontsize=14)
    
    # Save the figure
    unit_boxplot_file = output_dir / "unit_position_methylation_boxplot.png"
    plt.savefig(unit_boxplot_file, dpi=300, bbox_inches="tight")
    plt.close()
    visualizations['unit_position_boxplot'] = str(unit_boxplot_file)
    
    # Calculate and save unit methylation statistics
    unit_stats = {}
    for count, units in unit_methylation_by_count.items():
        unit_stats[count] = {}
        for unit_num, values in units.items():
            if values:
                unit_stats[count][unit_num] = {
                    'mean': float(np.mean(values)),
                    'median': float(np.median(values)),
                    'std': float(np.std(values)),
                    'count': len(values)
                }
    
    # Save statistics
    with open(stats_dir / "unit_methylation_stats.json", 'w') as f:
        json.dump(unit_stats, f, indent=2)
    
    logger.info(f"Unit methylation analysis complete. Results saved to {stats_dir}")
    
    return {
        'status': 'success',
        'unit_stats': unit_stats,
        'visualizations': visualizations
    }

def main():
    parser = argparse.ArgumentParser(description="Analyze upstream coding region methylation patterns")
    parser.add_argument("--data", "-d", required=True, help="Path to r-repeat methylation data pickle file")
    parser.add_argument("--output", "-o", help="Output directory for analysis results")
    parser.add_argument("--threshold", "-t", type=float, default=0.3, 
                      help="Threshold for low upstream methylation (default: 0.3)")
    
    args = parser.parse_args()
    
    analyze_upstream_coding_methylation(args.data, args.output, args.threshold)

if __name__ == "__main__":
    main()
