"""
Batch analysis script for downstream IGS region methylation in relation to r-repeats.

This script runs the downstream IGS methylation analysis across multiple samples
that already have r-repeat methylation data available.
"""

import argparse
import logging
import glob
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys

# Direct import from r_repeat package
from r_repeat.downstream_analysis.analyze_downstream_igs import analyze_downstream_igs_methylation

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler("logs/batch_downstream_analysis.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def batch_analyze_downstream_methylation(input_dir, output_dir=None, threshold=0.3):
    """
    Run downstream IGS methylation analysis on multiple samples.
    
    Parameters:
        input_dir: Directory containing r-repeat results with methylation pickle files
        output_dir: Parent directory for output files (default: same as input)
        threshold: Low methylation threshold (default: 0.3)
    
    Returns:
        dict: Dictionary with results for each sample
    """
    input_dir = Path(input_dir)
    
    if output_dir is None:
        output_dir = input_dir
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a batch data directory for cross-sample comparison
    batch_data_dir = output_dir / "batch_data"
    batch_data_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all r-repeat methylation pickle files by searching for data directories
    methylation_files = []
    
    # First check in standard data directory structure
    for data_dir in input_dir.glob("*/data"):
        meth_file = data_dir / "r_repeat_methylation.pkl"
        if meth_file.exists():
            methylation_files.append(meth_file)
    
    # Also look for files recursively in other locations (for backward compatibility)
    if not methylation_files:
        methylation_files = list(input_dir.glob("**/r_repeat_methylation.pkl"))
    
    if not methylation_files:
        logger.error(f"No r-repeat methylation files found in {input_dir}")
        return {"status": "error", "message": "No methylation files found"}
    
    # Process each methylation file
    results = {}
    batch_stats = {
        "sample_count": len(methylation_files),
        "low_downstream_threshold": threshold,
        "samples": {},
        "methylation_correlations": [],
        "downstream_stats": {
            "mean_values": [],
            "relation_to_r_repeat": {
                "mean_differences": [],
                "correlation_values": []
            },
            "low_downstream_counts": []
        },
        "unit_counts": {}
    }
    
    for meth_file in methylation_files:
        meth_path = Path(meth_file)
        
        # Determine sample name from directory structure
        if meth_path.parent.name == "data":
            sample_name = meth_path.parent.parent.name
        else:
            sample_name = meth_path.parent.name
        
        logger.info(f"Processing sample: {sample_name}")
        
        try:
            # Create organized output directory structure
            sample_output_dir = output_dir / sample_name / "downstream_analysis"
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            
            # Create data subdirectory for serialized results
            sample_data_dir = sample_output_dir / "data"
            sample_data_dir.mkdir(parents=True, exist_ok=True)
            
            # Run downstream IGS analysis
            sample_results = analyze_downstream_igs_methylation(
                meth_file, 
                output_dir=sample_output_dir,
                low_downstream_threshold=threshold
            )
            
            status = "success" if sample_results and sample_results.get('status') == 'success' else "failed"
            
            # Record basic results
            results[sample_name] = {
                "status": status,
                "output_dir": str(sample_output_dir),
                "data_dir": str(sample_data_dir),
                "visualizations": list(sample_results['visualizations'].values()) if sample_results and 'visualizations' in sample_results else []
            }
            
            # Save detailed sample-specific data for cross-sample comparison
            if status == "success":
                # Extract key statistics
                pearson_r = sample_results['pearson_correlation']['r']
                pearson_p = sample_results['pearson_correlation']['p_value']
                downstream_avg = sample_results['downstream_avg_methylation']
                r_repeat_avg = sample_results['r_repeat_avg_methylation']
                mean_diff = sample_results['mean_difference']
                region_count = sample_results['region_count']
                
                # Record in results
                results[sample_name].update({
                    "correlation": {
                        "pearson_r": float(pearson_r),
                        "pearson_p": float(pearson_p)
                    },
                    "methylation_levels": {
                        "downstream_avg": float(downstream_avg),
                        "r_repeat_avg": float(r_repeat_avg),
                        "mean_diff": float(mean_diff)
                    },
                    "region_count": region_count
                })
                
                # Add to batch statistics
                batch_stats["samples"][sample_name] = {
                    "status": "success",
                    "region_count": region_count,
                    "correlation": {
                        "pearson_r": float(pearson_r),
                        "pearson_p": float(pearson_p)
                    }
                }
                
                batch_stats["methylation_correlations"].append({
                    "sample": sample_name,
                    "pearson_r": float(pearson_r),
                    "pearson_p": float(pearson_p),
                    "significant": pearson_p < 0.05
                })
                
                batch_stats["downstream_stats"]["mean_values"].append(float(downstream_avg))
                batch_stats["downstream_stats"]["relation_to_r_repeat"]["mean_differences"].append(float(mean_diff))
                batch_stats["downstream_stats"]["relation_to_r_repeat"]["correlation_values"].append(float(pearson_r))
                
                # Try to load detailed stats if available
                try:
                    comparison_data_file = sample_output_dir / "stats" / "downstream_comparison_data.json"
                    if comparison_data_file.exists():
                        with open(comparison_data_file, 'r') as f:
                            comparison_data = json.load(f)
                            
                            # Record low downstream count
                            batch_stats["downstream_stats"]["low_downstream_counts"].append({
                                "sample": sample_name,
                                "low_downstream_count": comparison_data.get("low_downstream_count", 0),
                                "region_count": comparison_data.get("region_count", 0)
                            })
                            
                            # Aggregate unit count distribution
                            unit_distribution = comparison_data.get("unit_count_distribution", {})
                            for unit_count, count in unit_distribution.items():
                                if unit_count not in batch_stats["unit_counts"]:
                                    batch_stats["unit_counts"][unit_count] = {}
                                
                                batch_stats["unit_counts"][unit_count][sample_name] = count
                except Exception as e:
                    logger.warning(f"Could not load detailed comparison data for {sample_name}: {str(e)}")
                
            else:
                batch_stats["samples"][sample_name] = {
                    "status": "failed",
                    "error": sample_results.get('message', "Unknown error")
                }
            
        except Exception as e:
            logger.error(f"Error processing {sample_name}: {str(e)}", exc_info=True)
            results[sample_name] = {"status": "failed", "error": str(e)}
            batch_stats["samples"][sample_name] = {"status": "failed", "error": str(e)}
    
    # Calculate aggregate statistics
    if batch_stats["downstream_stats"]["mean_values"]:
        batch_stats["downstream_stats"]["aggregate"] = {
            "mean": np.mean(batch_stats["downstream_stats"]["mean_values"]),
            "median": np.median(batch_stats["downstream_stats"]["mean_values"]),
            "std_dev": np.std(batch_stats["downstream_stats"]["mean_values"]),
            "min": min(batch_stats["downstream_stats"]["mean_values"]),
            "max": max(batch_stats["downstream_stats"]["mean_values"])
        }
    
    if batch_stats["downstream_stats"]["relation_to_r_repeat"]["correlation_values"]:
        batch_stats["downstream_stats"]["relation_to_r_repeat"]["aggregate"] = {
            "mean_correlation": np.mean(batch_stats["downstream_stats"]["relation_to_r_repeat"]["correlation_values"]),
            "mean_difference": np.mean(batch_stats["downstream_stats"]["relation_to_r_repeat"]["mean_differences"])
        }
    
    # Generate cross-sample visualizations
    visualization_files = generate_cross_sample_visualizations(batch_stats, batch_data_dir, results)
    batch_stats["visualizations"] = visualization_files
    
    # Save the overall batch summary
    summary_file = batch_data_dir / "batch_downstream_analysis_summary.json"
    try:
        with open(summary_file, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Saved batch summary to {summary_file}")
    except Exception as e:
        logger.error(f"Error saving summary: {str(e)}", exc_info=True)
    
    # Save the detailed batch statistics
    stats_file = batch_data_dir / "batch_downstream_analysis_stats.json"
    try:
        with open(stats_file, 'w') as f:
            json.dump(batch_stats, f, indent=2)
        logger.info(f"Saved batch statistics to {stats_file}")
    except Exception as e:
        logger.error(f"Error saving batch statistics: {str(e)}", exc_info=True)
    
    logger.info(f"Batch analysis complete. Processed {len(methylation_files)} samples.")
    return results

def generate_cross_sample_visualizations(batch_stats, output_dir, results):
    """
    Generate visualizations comparing downstream methylation patterns across samples.
    
    Parameters:
        batch_stats: Dictionary with batch statistics
        output_dir: Directory for output files
        results: Dictionary with analysis results for each sample
        
    Returns:
        dict: Dictionary mapping visualization types to file paths
    """
    output_dir = Path(output_dir)
    visualizations_dir = output_dir / "visualizations"
    visualizations_dir.mkdir(parents=True, exist_ok=True)
    
    visualizations = {}
    
    # 1. Bar plot of correlation coefficients across samples
    if batch_stats["methylation_correlations"]:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Sort by correlation strength
        sorted_data = sorted(batch_stats["methylation_correlations"], 
                           key=lambda x: x["pearson_r"])
        
        samples = [item["sample"] for item in sorted_data]
        correlations = [item["pearson_r"] for item in sorted_data]
        significant = [item["significant"] for item in sorted_data]
        
        # Create colored bars based on significance
        colors = ['darkred' if sig else 'gray' for sig in significant]
        bars = ax.bar(samples, correlations, color=colors, alpha=0.8)
        
        # Add horizontal line at 0
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        
        # Add average correlation line
        if correlations:
            avg_corr = np.mean(correlations)
            ax.axhline(y=avg_corr, color='blue', linestyle='--', 
                      label=f'Mean correlation: {avg_corr:.3f}')
        
        ax.set_title("Correlation Between Downstream IGS and r-repeat Methylation Across Samples", 
                    fontsize=14)
        ax.set_xlabel("Sample", fontsize=12)
        ax.set_ylabel("Pearson Correlation Coefficient", fontsize=12)
        ax.set_ylim(-1, 1)
        ax.grid(axis='y', alpha=0.3)
        plt.xticks(rotation=90)
        ax.legend()
        
        # Save figure
        corr_file = visualizations_dir / "cross_sample_correlation.png"
        plt.savefig(corr_file, dpi=300, bbox_inches="tight")
        plt.close()
        visualizations["correlation_comparison"] = str(corr_file)
    
    # 2. Scatter plot of downstream vs r-repeat average methylation by sample
    if ("samples" in batch_stats and results and
        any("methylation_levels" in sample_data 
            for sample_name, sample_data in results.items() if "methylation_levels" in sample_data)):
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        samples = []
        downstream_values = []
        r_repeat_values = []
        
        for sample_name, sample_data in results.items():
            if "methylation_levels" in sample_data:
                samples.append(sample_name)
                downstream_values.append(sample_data["methylation_levels"]["downstream_avg"])
                r_repeat_values.append(sample_data["methylation_levels"]["r_repeat_avg"])
        
        if samples:
            # Create scatter plot
            ax.scatter(downstream_values, r_repeat_values, alpha=0.7)
            
            # Add labels to points
            for i, sample in enumerate(samples):
                ax.annotate(sample, (downstream_values[i], r_repeat_values[i]), 
                           fontsize=8, alpha=0.7)
            
            # Add linear regression if we have enough points
            if len(samples) > 2:
                from scipy.stats import linregress
                slope, intercept, r_value, p_value, _ = linregress(downstream_values, r_repeat_values)
                
                x_line = np.array([min(downstream_values), max(downstream_values)])
                y_line = slope * x_line + intercept
                ax.plot(x_line, y_line, 'r--', 
                       label=f'y={slope:.2f}x+{intercept:.2f} (r={r_value:.2f}, p={p_value:.4f})')
            
            ax.set_title("Downstream IGS vs r-repeat Methylation Across Samples", fontsize=14)
            ax.set_xlabel("Average Downstream Methylation", fontsize=12)
            ax.set_ylabel("Average r-repeat Methylation", fontsize=12)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.grid(True, alpha=0.3)
            if len(samples) > 2:
                ax.legend()
            
            # Save figure
            scatter_file = visualizations_dir / "downstream_vs_r_repeat_scatter.png"
            plt.savefig(scatter_file, dpi=300, bbox_inches="tight")
            plt.close()
            visualizations["downstream_vs_r_repeat_scatter"] = str(scatter_file)
    
    # 3. Heatmap of low downstream counts across samples and unit counts (if data available)
    # (Similar implementation as in the upstream analysis)
    
    logger.info(f"Generated {len(visualizations)} cross-sample visualizations in {visualizations_dir}")
    return visualizations

def main():
    parser = argparse.ArgumentParser(description="Batch analyze downstream IGS region methylation")
    parser.add_argument("--input", "-i", required=True, 
                        help="Directory containing r-repeat results with methylation pickle files")
    parser.add_argument("--output", "-o", 
                        help="Parent directory for output files")
    parser.add_argument("--threshold", "-t", type=float, default=0.3,
                        help="Threshold for low downstream methylation (default: 0.3)")
    
    args = parser.parse_args()
    
    batch_analyze_downstream_methylation(args.input, args.output, args.threshold)

if __name__ == "__main__":
    main()
