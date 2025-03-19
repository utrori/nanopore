"""
Analysis module for downstream intergenic spacer (IGS) region methylation patterns in relation to r-repeat regions.

This module is designed to be used as part of the r-repeat methylation analysis pipeline
but can also be run independently for more detailed downstream IGS region analysis.
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import argparse
import logging
import json

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

def analyze_downstream_igs_methylation(data_file, output_dir=None, low_downstream_threshold=0.3):
    """
    Analyze the relationship between downstream IGS region methylation and r-repeat methylation.
    
    Parameters:
        data_file: Path to a pickle file containing r-repeat methylation data with downstream information
        output_dir: Directory for output files
        low_downstream_threshold: Threshold to define low downstream methylation (default: 0.3)
        
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
        output_dir = parent_dir.parent / "downstream_analysis"
    
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
    downstream_meth = []
    r_repeat_meth = []
    region_count = 0
    
    # Lists to store regions with low downstream methylation
    low_downstream_regions = []
    
    if has_multiple_regions:
        # Handle multiple regions per read
        for read_id, regions in r_repeat_data.items():
            for region in regions:
                if ('downstream_igs_methylation' in region and 
                    region['downstream_igs_methylation']['available']):
                    
                    downstream_avg = region['downstream_igs_methylation']['avg_methylation'] / 255
                    r_repeat_avg = region['avg_methylation'] / 255
                    
                    if downstream_avg >= 0:
                        downstream_meth.append(downstream_avg)
                        r_repeat_meth.append(r_repeat_avg)
                        region_count += 1
                        
                        # Store regions with low downstream methylation
                        if downstream_avg < low_downstream_threshold:
                            low_downstream_regions.append({
                                'read_id': read_id,
                                'region_data': region,
                                'downstream_methylation': downstream_avg,
                                'r_repeat_methylation': r_repeat_avg
                            })
    else:
        # Handle single region per read
        for read_id, data in r_repeat_data.items():
            if ('downstream_igs_methylation' in data and 
                data['downstream_igs_methylation']['available']):
                
                downstream_avg = data['downstream_igs_methylation']['avg_methylation'] / 255
                r_repeat_avg = data['avg_methylation'] / 255
                
                if downstream_avg >= 0:
                    downstream_meth.append(downstream_avg)
                    r_repeat_meth.append(r_repeat_avg)
                    region_count += 1
                    
                    # Store regions with low downstream methylation
                    if downstream_avg < low_downstream_threshold:
                        low_downstream_regions.append({
                            'read_id': read_id,
                            'region_data': data,
                            'downstream_methylation': downstream_avg,
                            'r_repeat_methylation': r_repeat_avg
                        })
    
    logger.info(f"Analyzed {region_count} regions with available downstream IGS data")
    logger.info(f"Found {len(low_downstream_regions)} regions with low downstream methylation (<{low_downstream_threshold})")
    
    if len(downstream_meth) < 2:
        logger.warning("Not enough data points with downstream IGS information for analysis")
        return {
            'status': 'error',
            'message': 'Not enough data points with downstream IGS information',
            'visualizations': {}
        }
    
    # Calculate basic statistics
    downstream_avg = np.mean(downstream_meth)
    r_repeat_avg = np.mean(r_repeat_meth)
    meth_diffs = [r - d for r, d in zip(r_repeat_meth, downstream_meth)]
    avg_diff = np.mean(meth_diffs)
    
    # Calculate correlation
    from scipy.stats import pearsonr, spearmanr
    pearson_r, pearson_p = pearsonr(downstream_meth, r_repeat_meth)
    spearman_r, spearman_p = spearmanr(downstream_meth, r_repeat_meth)
    
    # Perform linear regression
    from scipy.stats import linregress
    slope, intercept, r_value, p_value, std_err = linregress(downstream_meth, r_repeat_meth)
    
    # Plot 1: Scatter plot with regression line
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(downstream_meth, r_repeat_meth, alpha=0.6)
    
    # Add regression line
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
    ax.legend()
    
    scatter_file = output_dir / "downstream_vs_r_repeat_scatter.png"
    plt.savefig(scatter_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Plot 2: Distribution of methylation differences
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(meth_diffs, bins=30, kde=True, ax=ax)
    
    ax.axvline(x=0, color='r', linestyle='--', 
              label='No difference between downstream and r-repeat')
    ax.axvline(x=avg_diff, color='g', linestyle='-', 
              label=f'Mean difference: {avg_diff:.3f}')
    
    ax.set_title("Distribution of Methylation Differences (r-repeat - Downstream)", fontsize=14)
    ax.set_xlabel("Methylation Difference", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.legend()
    
    diff_file = output_dir / "methylation_difference_distribution.png"
    plt.savefig(diff_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Plot 3: Heatmap of methylation status (binned)
    # Convert continuous methylation levels to binned categories for visualization
    df = pd.DataFrame({
        'Downstream': downstream_meth,
        'r-repeat': r_repeat_meth
    })
    
    # Create bins for methylation levels
    bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    labels = ['0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0']
    
    df['Downstream_binned'] = pd.cut(df['Downstream'], bins=bins, labels=labels)
    df['r-repeat_binned'] = pd.cut(df['r-repeat'], bins=bins, labels=labels)
    
    # Create contingency table
    ctab = pd.crosstab(df['Downstream_binned'], df['r-repeat_binned'], normalize='all') * 100
    
    # Plot heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(ctab, annot=True, cmap="YlGnBu", fmt=".1f", cbar_kws={'label': 'Percentage (%)'})
    
    ax.set_title("Relationship Between Downstream IGS and r-repeat Methylation Levels", fontsize=14)
    ax.set_xlabel("r-repeat Region Methylation", fontsize=12)
    ax.set_ylabel("Downstream IGS Region Methylation", fontsize=12)
    
    heatmap_file = output_dir / "methylation_relationship_heatmap.png"
    plt.savefig(heatmap_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Create summary statistics
    summary = {
        'status': 'success',
        'downstream_avg_methylation': downstream_avg,
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
    
    # Create comprehensive serialized data for cross-sample comparison
    comparison_data = {
        "sample_name": sample_name,
        "region_count": region_count,
        "low_downstream_count": len(low_downstream_regions),
        "low_downstream_threshold": low_downstream_threshold,
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
            "downstream": {
                "mean": float(downstream_avg),
                "quartiles": {
                    "q1": float(np.percentile(downstream_meth, 25)),
                    "median": float(np.median(downstream_meth)),
                    "q3": float(np.percentile(downstream_meth, 75))
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
    for region_info in low_downstream_regions:
        region_data = region_info['region_data']
        if 'coordinates' in region_data and 'unit_count' in region_data['coordinates']:
            unit_count = region_data['coordinates']['unit_count']
            if unit_count not in unit_count_data:
                unit_count_data[unit_count] = 0
            unit_count_data[unit_count] += 1
    
    comparison_data["unit_count_distribution"] = {str(k): v for k, v in unit_count_data.items()}
    
    # Save comprehensive comparison data
    comparison_file = stats_dir / "downstream_comparison_data.json"
    with open(comparison_file, 'w') as f:
        json.dump(comparison_data, f, indent=2)
    
    # Save summary to file
    summary_file = output_dir / "downstream_analysis_summary.json"
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
    logger.info(f"Mean difference (r-repeat - downstream): {avg_diff:.3f}")
    logger.info(f"Cross-sample comparison data saved to {comparison_file}")
    
    return summary

def main():
    parser = argparse.ArgumentParser(description="Analyze downstream IGS region methylation patterns")
    parser.add_argument("--data", "-d", required=True, help="Path to r-repeat methylation data pickle file")
    parser.add_argument("--output", "-o", help="Output directory for analysis results")
    parser.add_argument("--threshold", "-t", type=float, default=0.3, 
                      help="Threshold for low downstream methylation (default: 0.3)")
    
    args = parser.parse_args()
    
    analyze_downstream_igs_methylation(args.data, args.output, args.threshold)

if __name__ == "__main__":
    main()
