"""
R-repeat workflow management script.

This script coordinates the complete r-repeat analysis pipeline by integrating:
1. r_repeat_detection.py - Detecting r-repeat regions
2. r_repeat_methylation_analysis.py - Analyzing methylation patterns
3. r_repeat_unit_analysis_new.py - Analyzing internal unit structure (NEW)
4. r_repeat_multi_analysis.py - Analyzing multi-region reads (when applicable)

It handles the extraction of reads to FASTQ format when necessary and manages
the flow of data between analysis steps.
"""

import os
import subprocess
import json
import numpy as np
import pickle
import argparse
import logging
from pathlib import Path
import shutil
import copy  # Add this import at the top of the file

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler("logs/r_repeat_workflow.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Import the other scripts as modules (no need to re-implement them)
import r_repeat_detection
import r_repeat_methylation_analysis
import r_repeat_unit_analysis_new  # Use the new unit analysis script
import r_repeat_multi_analysis
# Remove direct import of analyze_upstream_coding as it will be called via r_repeat_methylation_analysis

def extract_reads_to_fastq(bam_path, output_fastq=None, temp_dir=None):
    """
    Extract reads from a BAM file to FASTQ format.
    
    Parameters:
        bam_path: Path to BAM file
        output_fastq: Optional path for output FASTQ file
        temp_dir: Optional directory for temporary files
    
    Returns:
        str: Path to the output FASTQ file
    """
    # Convert input path to Path object if it's a string
    bam_path = Path(bam_path)
    
    # Create a temporary file if no output path provided
    if output_fastq is None:
        if temp_dir is None:
            temp_dir = "temp_files"
        
        temp_dir = Path(temp_dir)
        temp_dir.mkdir(parents=True, exist_ok=True)
        output_fastq = temp_dir / f"{bam_path.stem}.fastq"
    
    # Ensure output_fastq is a Path object
    output_fastq = Path(output_fastq)
    
    # Skip extraction if file already exists and has content
    if output_fastq.exists() and output_fastq.stat().st_size > 0:
        logger.info(f"Using existing FASTQ file: {output_fastq}")
        return str(output_fastq)
    
    # Convert BAM to FASTQ using samtools
    logger.info(f"Extracting reads from {bam_path} to {output_fastq}")
    
    cmd = ["samtools", "fastq", "-@", "4", str(bam_path)]
    with open(str(output_fastq), 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)
    
    return str(output_fastq)

def run_complete_analysis(bam_path, output_dir=None, keep_temp_files=False, reuse_existing=False, allow_multiple=True):
    """
    Run complete r-repeat analysis workflow on a single BAM file.
    
    This function coordinates the full analysis pipeline:
    1. Extract reads to FASTQ format
    2. Detect r-repeat regions
    3. Analyze unit structure with the new r_repeat_unit_analysis_new script
    4. Analyze multi-region patterns (if allow_multiple=True)
    5. Analyze methylation patterns (including upstream coding and downstream IGS analysis)
    
    Parameters:
        bam_path: Path to BAM file with nanopore reads (Dorado-called with methylation)
        output_dir: Directory for all output files
        keep_temp_files: Whether to keep temporary files after analysis
        reuse_existing: Whether to reuse existing intermediate files
        allow_multiple: Whether to allow multiple r-repeat regions per read
        
    Returns:
        dict: Results summary containing paths to all output files
    """
    bam_path = Path(bam_path)
    
    # Set up output directory with proper sample name
    sample_name = bam_path.stem
    if output_dir is None:
        output_dir = Path("r_repeat_results") / sample_name
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up temporary directory inside the sample directory
    temp_dir = output_dir / "temp"
    temp_dir.mkdir(exist_ok=True)
    
    # Create dedicated output subdirectories for different types of analyses
    unit_dir = output_dir / "unit_analysis"
    unit_dir.mkdir(exist_ok=True)
    
    methylation_dir = output_dir / "methylation"
    methylation_dir.mkdir(exist_ok=True)
    
    # Upstream and downstream analysis directories are part of methylation analysis
    upstream_dir = methylation_dir / "upstream_analysis"
    upstream_dir.mkdir(exist_ok=True)
    
    downstream_dir = methylation_dir / "downstream_analysis"
    downstream_dir.mkdir(exist_ok=True)
    
    multi_region_dir = output_dir / "multi_region_analysis"
    multi_region_dir.mkdir(exist_ok=True)
    
    # Create directory for serialized data
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)
    
    try:
        # Step 1: Extract reads to FASTQ (will reuse existing file if available)
        fastq_path = extract_reads_to_fastq(bam_path, temp_dir=temp_dir)
        logger.info(f"Using FASTQ: {fastq_path}")
        
        # Step 2: Detect r-repeat regions using reuse_existing flag
        logger.info("Step 1/5: Detecting r-repeat regions...")
        
        # Run detection pipeline with reuse flag
        r_repeat_regions, r_repeat_fasta = r_repeat_detection.detect_r_repeats(
            fastq_path, 
            output_dir=output_dir,
            reuse_existing=reuse_existing,
            allow_multiple=allow_multiple,
            temp_dir=temp_dir  # Pass temp directory for alignment files
        )
        
        # Check if r-repeat regions were found
        if not r_repeat_regions:
            logger.warning("No r-repeat regions found. Analysis stopping.")
            return {"status": "failed", "error": "No r-repeat regions found"}

        # Step 3: Analyze unit structure
        logger.info("Step 2/5: Analyzing r-repeat unit structure...")
        
        # Create motif reference
        motif_ref = r_repeat_unit_analysis_new.create_motif_reference()
        
        # Align r-repeat sequences to motif
        alignment_bam = r_repeat_unit_analysis_new.align_to_motif(
            fastq_path=fastq_path, 
            output_dir=temp_dir,  # Store alignments in temp directory
            reuse_existing=reuse_existing
        )
        
        # Process motif alignments to identify units
        read_alignments = r_repeat_unit_analysis_new.process_motif_alignments(alignment_bam)
        
        # Add correspondence between units and regions, and get distributions
        unit_count_distribution, unit_lengths, region_lengths, read_ids_with_problematic_units, unit_analysis_data = r_repeat_unit_analysis_new.correspond_units_to_regions(
            read_alignments, r_repeat_regions)
        
        # Generate visualizations in dedicated unit analysis directory
        unit_count_viz = r_repeat_unit_analysis_new.generate_unit_count_visualization(
            unit_count_distribution, output_dir=unit_dir
        )
        
        unit_length_viz = r_repeat_unit_analysis_new.generate_unit_length_visualization(
            unit_lengths, output_dir=unit_dir
        )
        
        region_length_viz = r_repeat_unit_analysis_new.generate_region_length_visualization(
            region_lengths, output_dir=unit_dir
        )
        
        # Serialize unit analysis results
        unit_data = {
            'unit_count_distribution': unit_count_distribution,
            'unit_lengths': unit_lengths,
            'region_lengths': region_lengths,
            'read_ids_with_problematic_units': list(read_ids_with_problematic_units),
            'sample_name': sample_name
        }
        
        # Add the comprehensive unit analysis data to unit_data
        unit_data.update(unit_analysis_data)
        
        with open(data_dir / "unit_analysis.json", 'w') as f:
            json.dump(unit_data, f, indent=2)
        
        # Step 4: Multi-region analysis (if allowed)
        multi_region_results = None
        if allow_multiple:
            logger.info("Step 3/5: Analyzing multi-region patterns...")
            
            # Create a serializable version of r_repeat_regions that includes unit information
            # Keep regions with problematic units in JSON, but mark them with flags
            serializable_regions = {}
            for read_id, regions in r_repeat_regions.items():
                if isinstance(regions, list):
                    # Multiple regions per read
                    serializable_regions[read_id] = []
                    for region in regions:
                        region_copy = {
                            'start': int(region['start']),
                            'end': int(region['end']),
                            'length': int(region['length']),
                            'is_reverse': bool(region['is_reverse']),
                            'unit_count': region.get('unit_count', 0),
                            'has_problematic_unit': region.get('has_problematic_unit', False),
                            'has_long_unit': region.get('has_long_unit', False),
                            'has_short_unit': region.get('has_short_unit', False)
                        }
                        
                        # Add unit lengths if available
                        if 'unit_lengths' in region:
                            region_copy['unit_lengths'] = [int(length) for length in region['unit_lengths']]
                            
                        serializable_regions[read_id].append(region_copy)
                else:
                    # Single region per read
                    region = regions
                    region_copy = {
                        'start': int(region['start']),
                        'end': int(region['end']),
                        'length': int(region['length']),
                        'is_reverse': bool(region['is_reverse']),
                        'unit_count': region.get('unit_count', 0),
                        'has_problematic_unit': region.get('has_problematic_unit', False),
                        'has_long_unit': region.get('has_long_unit', False),
                        'has_short_unit': region.get('has_short_unit', False)
                    }
                    
                    # Add unit lengths if available
                    if 'unit_lengths' in region:
                        region_copy['unit_lengths'] = [int(length) for length in region['unit_lengths']]
                        
                    serializable_regions[read_id] = region_copy
            
            # Save the enhanced regions data in the data directory
            regions_file = data_dir / "r_repeat_regions.json"
            with open(regions_file, 'w') as f:
                json.dump(serializable_regions, f, indent=2)
            
            # Run multi-region analysis
            multi_regions, multi_viz = r_repeat_multi_analysis.analyze_multi_regions(
                regions_file,
                output_dir=multi_region_dir
            )
            
            # Store results
            multi_region_results = {
                'multi_regions': multi_regions,
                'visualizations': multi_viz
            }
            
            # Save serialized multi-region analysis data
            with open(data_dir / "multi_region_analysis.json", 'w') as f:
                multi_data = {
                    'sample_name': sample_name,
                    'read_count': len(multi_regions),
                    'region_count': sum(len(regions) for regions in multi_regions.values()),
                    'multi_region_stats': r_repeat_multi_analysis.get_multi_region_statistics(multi_regions)
                }
                json.dump(multi_data, f, indent=2)
        else:
            # Save regions to file for single-region mode
            regions_file = data_dir / "r_repeat_regions.json"
            if not regions_file.exists():
                with open(regions_file, 'w') as f:
                    # Convert r-repeat regions to a serializable format with unit info
                    # Keep all regions but mark those with problematic units
                    serializable_regions = {}
                    for read_id, region in r_repeat_regions.items():
                        region_copy = {
                            'start': int(region['start']),
                            'end': int(region['end']),
                            'length': int(region['length']),
                            'is_reverse': bool(region['is_reverse']),
                            'unit_count': region.get('unit_count', 0),
                            'has_problematic_unit': region.get('has_problematic_unit', False),
                            'has_long_unit': region.get('has_long_unit', False),
                            'has_short_unit': region.get('has_short_unit', False)
                        }
                        
                        # Add unit lengths if available
                        if 'unit_lengths' in region:
                            region_copy['unit_lengths'] = [int(length) for length in region['unit_lengths']]
                            
                        serializable_regions[read_id] = region_copy
                    
                    json.dump(serializable_regions, f, indent=2)
        
        # Step 5: Analyze methylation patterns (now includes upstream coding and downstream IGS analysis)
        logger.info("Step 4/4: Analyzing methylation patterns (including upstream and downstream regions)...")
        
        # Create a clean copy of r_repeat_regions without any unpicklable objects
        r_repeat_regions_clean = {}
        for read_id, regions in r_repeat_regions.items():
            if isinstance(regions, list):
                # Multiple regions per read
                valid_regions = []
                for region in regions:
                    region_copy = region.copy()
                    if 'units' in region_copy:
                        del region_copy['units']
                    valid_regions.append(region_copy)
                r_repeat_regions_clean[read_id] = valid_regions
            else:
                # Single region per read
                region = regions
                region_copy = region.copy()
                if 'units' in region_copy:
                    del region_copy['units']
                r_repeat_regions_clean[read_id] = region_copy
        
        # Run methylation analysis with the cleaned data structure
        # This will now internally handle both upstream and downstream region analysis
        r_repeat_data, methylation_viz, upstream_results, downstream_results = r_repeat_methylation_analysis.analyze_methylation(
            r_repeat_regions_clean,
            bam_path,
            output_dir=methylation_dir,
            analyze_upstream=True,  # Flag to enable upstream analysis
            analyze_downstream=True  # Flag to enable downstream analysis
        )
        
        # Extract upstream visualization paths from the results
        upstream_viz = []
        if isinstance(upstream_results, dict) and 'visualizations' in upstream_results:
            upstream_viz = list(upstream_results['visualizations'].values())
            
        # Extract downstream visualization paths from the results
        downstream_viz = []
        if isinstance(downstream_results, dict) and 'visualizations' in downstream_results:
            downstream_viz = list(downstream_results['visualizations'].values())
        
        # Compile comprehensive results summary
        results = {
            "status": "success",
            "sample_name": sample_name,
            "output_directory": str(output_dir),
            "directories": {
                "temp": str(temp_dir),
                "unit_analysis": str(unit_dir),
                "methylation": str(methylation_dir),
                "upstream_analysis": str(upstream_dir),
                "downstream_analysis": str(downstream_dir),
                "multi_region_analysis": str(multi_region_dir) if allow_multiple else None,
                "data": str(data_dir)
            },
            "files": {
                "r_repeat_regions_json": str(data_dir / "r_repeat_regions.json"),
                "r_repeat_sequences_fasta": r_repeat_fasta,
                "r_repeat_methylation": str(data_dir / "r_repeat_methylation.pkl"),
                "unit_analysis_json": str(data_dir / "unit_analysis.json"),
                "motif_alignment_bam": alignment_bam
            },
            "visualizations": {
                "unit_count": unit_count_viz,
                "unit_length": unit_length_viz,
                "region_length": region_length_viz,
                "methylation": methylation_viz,
                "upstream": upstream_viz,
                "downstream": downstream_viz,  # Add the downstream viz from the methylation analysis
                "multi_region": multi_region_results["visualizations"] if multi_region_results else []
            },
            "allow_multiple": allow_multiple,
            "unit_count_distribution": unit_count_distribution,
            "excluded_read_count": len(read_ids_with_problematic_units),
            "unit_lengths_summary": {
                "count": len(unit_lengths),
                "min": min(unit_lengths) if unit_lengths else None,
                "mean": float(np.mean(unit_lengths)) if unit_lengths else None, 
                "median": float(np.median(unit_lengths)) if unit_lengths else None,
                "max": max(unit_lengths) if unit_lengths else None
            },
            "unit_filter_parameters": {
                "min_unit_length": 300,
                "max_unit_length": 1000
            },
            "region_lengths_summary": {
                "count": len(region_lengths),
                "min": min(region_lengths) if region_lengths else None,
                "mean": float(np.mean(region_lengths)) if region_lengths else None, 
                "median": float(np.median(region_lengths)) if region_lengths else None,
                "max": max(region_lengths) if region_lengths else None
            }
        }
        
        if allow_multiple:
            # Count total regions and calculate statistics
            region_count = sum(len(regions) if isinstance(regions, list) else 1 
                              for regions in r_repeat_regions.values())
            multi_region_reads = sum(1 for regions in r_repeat_regions.values() 
                                    if isinstance(regions, list) and len(regions) > 1)
            
            # Count reads by unit count
            unit_counts = {}
            for read_id, regions in r_repeat_regions.items():
                if isinstance(regions, list):
                    for region in regions:
                        unit_count = region.get('unit_count', 0)
                        unit_counts[unit_count] = unit_counts.get(unit_count, 0) + 1
                else:
                    unit_count = regions.get('unit_count', 0)
                    unit_counts[unit_count] = unit_counts.get(unit_count, 0) + 1
            
            # Add statistics to results
            results.update({
                "r_repeat_count": region_count,
                "read_count": len(r_repeat_regions),
                "multi_region_reads": multi_region_reads,
                "unit_counts": unit_counts
            })
        else:
            results["r_repeat_count"] = len(r_repeat_regions)
            
            # Count reads by unit count for single-region mode
            unit_counts = {}
            for read_id, region in r_repeat_regions.items():
                unit_count = region.get('unit_count', 0)
                unit_counts[unit_count] = unit_counts.get(unit_count, 0) + 1
            
            results["unit_counts"] = unit_counts
        
        # Save summary to file
        summary_file = data_dir / "analysis_summary.json"
        with open(summary_file, 'w') as f:
            # Convert paths to strings for JSON serialization
            json_results = {
                k: (v if not isinstance(v, Path) else str(v))
                for k, v in results.items()
            }
            
            # Handle nested dictionaries and lists with Path objects
            if isinstance(json_results["directories"], dict):
                json_results["directories"] = {
                    k: (v if not isinstance(v, Path) else str(v))
                    for k, v in json_results["directories"].items()
                }
                
            if isinstance(json_results["files"], dict):
                json_results["files"] = {
                    k: (v if not isinstance(v, Path) else str(v))
                    for k, v in json_results["files"].items()
                }
                
            if isinstance(json_results["visualizations"], dict):
                json_results["visualizations"] = {
                    k: ([str(item) for item in v] if isinstance(v, list) else
                        (str(v) if isinstance(v, Path) else v))
                    for k, v in json_results["visualizations"].items()
                }
            
            json.dump(json_results, f, indent=2)
        
        # Create a symlink to the latest results for convenience
        latest_link = Path("r_repeat_results") / "latest"
        if latest_link.is_symlink() or latest_link.exists():
            latest_link.unlink()
        
        try:
            latest_link.symlink_to(output_dir)
        except:
            logger.warning("Could not create 'latest' symlink. This may be due to filesystem limitations.")
        
        logger.info(f"Analysis complete. Found {results.get('r_repeat_count', 0)} r-repeat regions.")
        logger.info(f"Results saved to {output_dir}")
        
        return results
    finally:
        # Clean up temporary files unless instructed to keep them
        if not keep_temp_files:
            try:
                if Path(fastq_path).exists():
                    Path(fastq_path).unlink()
                    logger.info(f"Removed temporary FASTQ file: {fastq_path}")
                
                # Remove other temporary files but keep the directory structure
                for temp_file in temp_dir.glob("*"):
                    if temp_file.is_file():
                        temp_file.unlink()
            except Exception as e:
                logger.warning(f"Error cleaning up temporary files: {str(e)}")

def run_batch_analysis(input_dir, output_parent_dir=None, pattern="*.bam", keep_temp_files=False, reuse_existing=False, allow_multiple=True):
    """
    Run r-repeat analysis on multiple BAM files in a directory.
    
    Parameters:
        input_dir: Directory containing BAM files
        output_parent_dir: Parent directory for all output directories
        pattern: File pattern to match BAM files
        keep_temp_files: Whether to keep temporary files after analysis
        reuse_existing: Whether to reuse existing intermediate files when possible
        allow_multiple: Whether to allow multiple r-repeat regions per read
        
    Returns:
        dict: Results summary for all samples
    """
    input_dir = Path(input_dir)
    
    if output_parent_dir is None:
        output_parent_dir = Path("r_repeat_results")
    
    output_parent_dir = Path(output_parent_dir)
    output_parent_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all BAM files
    bam_files = list(input_dir.glob(pattern))
    
    if not bam_files:
        logger.warning(f"No BAM files found in {input_dir} matching pattern '{pattern}'")
        return {"status": "failed", "error": "No BAM files found"}
    
    logger.info(f"Found {len(bam_files)} BAM files to process")
    
    # Process each BAM file
    results = {}
    
    # Add batch summary serialization for cross-sample comparison
    batch_summary = {
        "status": "success",
        "samples": {},
        "sample_count": len(bam_files),
        "r_repeat_total_count": 0,
        "unit_statistics": {
            "count_distribution": {},
            "length_statistics": {
                "mean": [],
                "median": [],
                "min": [],
                "max": []
            }
        },
        "upstream_methylation": {
            "correlation": [],
            "mean_difference": []
        }
    }
    
    for bam_file in bam_files:
        sample_name = bam_file.stem
        logger.info(f"Processing sample: {sample_name}")
        
        try:
            sample_output_dir = output_parent_dir / sample_name
            sample_results = run_complete_analysis(
                bam_file,
                output_dir=sample_output_dir,
                keep_temp_files=keep_temp_files,
                reuse_existing=reuse_existing,
                allow_multiple=allow_multiple
            )
            
            # Add sample results to the batch summary
            batch_summary["samples"][sample_name] = {
                "status": sample_results.get("status", "failed"),
                "r_repeat_count": sample_results.get("r_repeat_count", 0),
                "unit_counts": sample_results.get("unit_counts", {}),
                "unit_lengths_summary": sample_results.get("unit_lengths_summary", {})
            }
            
            # Collect statistics for cross-sample analysis
            if sample_results.get("status") == "success":
                # Add to total r-repeat count
                batch_summary["r_repeat_total_count"] += sample_results.get("r_repeat_count", 0)
                
                # Aggregate unit count distribution
                for count, freq in sample_results.get("unit_counts", {}).items():
                    count_str = str(count)  # Convert to string for JSON keys
                    if count_str not in batch_summary["unit_statistics"]["count_distribution"]:
                        batch_summary["unit_statistics"]["count_distribution"][count_str] = 0
                    batch_summary["unit_statistics"]["count_distribution"][count_str] += freq
                
                # Add unit length statistics
                if "unit_lengths_summary" in sample_results and sample_results["unit_lengths_summary"].get("count", 0) > 0:
                    batch_summary["unit_statistics"]["length_statistics"]["mean"].append(sample_results["unit_lengths_summary"].get("mean"))
                    batch_summary["unit_statistics"]["length_statistics"]["median"].append(sample_results["unit_lengths_summary"].get("median"))
                    batch_summary["unit_statistics"]["length_statistics"]["min"].append(sample_results["unit_lengths_summary"].get("min"))
                    batch_summary["unit_statistics"]["length_statistics"]["max"].append(sample_results["unit_lengths_summary"].get("max"))
                
                # Try to load upstream methylation summary and add to batch statistics
                try:
                    upstream_summary_file = sample_output_dir / "upstream_analysis" / "upstream_analysis_summary.json"
                    if upstream_summary_file.exists():
                        with open(upstream_summary_file, 'r') as f:
                            upstream_data = json.load(f)
                            
                            if "pearson_correlation" in upstream_data:
                                batch_summary["upstream_methylation"]["correlation"].append({
                                    "sample": sample_name,
                                    "pearson_r": upstream_data["pearson_correlation"]["r"],
                                    "p_value": upstream_data["pearson_correlation"]["p_value"]
                                })
                                
                            if "mean_difference" in upstream_data:
                                batch_summary["upstream_methylation"]["mean_difference"].append({
                                    "sample": sample_name,
                                    "value": upstream_data["mean_difference"]
                                })
                except Exception as e:
                    logger.warning(f"Error loading upstream analysis for {sample_name}: {str(e)}")
            
            results[sample_name] = sample_results
            
        except Exception as e:
            logger.error(f"Error processing {sample_name}: {str(e)}", exc_info=True)
            results[sample_name] = {"status": "failed", "error": str(e)}
            batch_summary["samples"][sample_name] = {"status": "failed", "error": str(e)}
    
    # Calculate aggregate statistics for the batch
    if batch_summary["unit_statistics"]["length_statistics"]["mean"]:
        batch_summary["unit_statistics"]["aggregate"] = {
            "mean_length": np.mean(batch_summary["unit_statistics"]["length_statistics"]["mean"]),
            "median_length": np.median(batch_summary["unit_statistics"]["length_statistics"]["median"]),
            "min_length": min(batch_summary["unit_statistics"]["length_statistics"]["min"]),
            "max_length": max(batch_summary["unit_statistics"]["length_statistics"]["max"])
        }
    
    # Save the detailed batch summary for cross-sample analysis
    batch_data_dir = output_parent_dir / "batch_data"
    batch_data_dir.mkdir(exist_ok=True)
    
    batch_summary_file = batch_data_dir / "batch_analysis_summary.json"
    with open(batch_summary_file, 'w') as f:
        json.dump(batch_summary, f, indent=2)
    
    # Save cross-sample comparison visualizations
    generate_cross_sample_visualizations(batch_summary, batch_data_dir)
    
    logger.info(f"Batch analysis complete. Processed {len(bam_files)} samples.")
    logger.info(f"Batch summary saved to {batch_summary_file}")
    
    return results

def generate_cross_sample_visualizations(batch_summary, output_dir):
    """
    Generate visualizations comparing results across different samples.
    
    Parameters:
        batch_summary: Dictionary with batch analysis summary data
        output_dir: Directory for output files
        
    Returns:
        list: List of paths to generated visualization files
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    output_files = []
    
    # Figure 1: Unit count distribution across samples
    if batch_summary["unit_statistics"]["count_distribution"]:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Convert string keys to int and sort
        unit_counts = sorted(int(k) for k in batch_summary["unit_statistics"]["count_distribution"].keys())
        frequencies = [batch_summary["unit_statistics"]["count_distribution"][str(count)] for count in unit_counts]
        
        bars = ax.bar(unit_counts, frequencies)
        
        # Add count labels on top of bars
        for bar, freq in zip(bars, frequencies):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{freq}', ha='center', va='bottom')
        
        ax.set_title("Distribution of r-repeat Units per Region Across All Samples", fontsize=14)
        ax.set_xlabel("Number of Units", fontsize=12)
        ax.set_ylabel("Total Count", fontsize=12)
        ax.set_xticks(unit_counts)
        ax.grid(axis='y', alpha=0.3)
        
        # Save figure
        output_file = output_dir / "cross_sample_unit_distribution.png"
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(output_file))
    
    # Figure 2: Upstream methylation correlation coefficients across samples
    if batch_summary["upstream_methylation"]["correlation"]:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Extract data for plotting
        samples = [item["sample"] for item in batch_summary["upstream_methylation"]["correlation"]]
        correlations = [item["pearson_r"] for item in batch_summary["upstream_methylation"]["correlation"]]
        p_values = [item["p_value"] for item in batch_summary["upstream_methylation"]["correlation"]]
        
        # Sort by correlation strength
        sorted_indices = np.argsort(correlations)
        samples = [samples[i] for i in sorted_indices]
        correlations = [correlations[i] for i in sorted_indices]
        p_values = [p_values[i] for i in sorted_indices]
        
        # Create bars with color based on statistical significance (p < 0.05)
        bars = ax.bar(samples, correlations, 
                     color=['red' if p < 0.05 else 'gray' for p in p_values])
        
        # Add horizontal line at 0
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        
        ax.set_title("Correlation Between Upstream Coding and r-repeat Methylation by Sample", fontsize=14)
        ax.set_xlabel("Sample", fontsize=12)
        ax.set_ylabel("Pearson Correlation Coefficient", fontsize=12)
        ax.set_ylim(-1, 1)
        ax.grid(axis='y', alpha=0.3)
        plt.xticks(rotation=90)
        
        # Save figure
        output_file = output_dir / "cross_sample_upstream_correlation.png"
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()
        output_files.append(str(output_file))
    
    logger.info(f"Generated {len(output_files)} cross-sample comparison visualizations")
    return output_files

def main():
    """Main entry point for command-line execution"""
    parser = argparse.ArgumentParser(
        description="Run complete r-repeat analysis workflow on nanopore data"
    )
    
    # Define input mode options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--bam", "-b", 
        help="Path to single BAM file with nanopore reads"
    )
    input_group.add_argument(
        "--dir", "-d", 
        help="Directory containing multiple BAM files to analyze"
    )
    
    # Additional options
    parser.add_argument(
        "--output", "-o",
        help="Output directory for analysis results"
    )
    parser.add_argument(
        "--pattern", "-p",
        default="*.bam",
        help="File pattern for batch mode (default: *.bam)"
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary files after analysis"
    )
    
    # Add reuse option
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="Reuse existing intermediate files when available (requires --keep-temp from previous run)"
    )
    
    # Add option to allow multiple r-repeat regions
    parser.add_argument(
        "--multiple",
        action="store_true",
        help="Allow detection of multiple r-repeat regions per read"
    )
    
    # Parse arguments and run analysis
    args = parser.parse_args()
    
    # Create log directory if it doesn't exist
    Path("logs").mkdir(exist_ok=True)
    
    if args.bam:
        # Single file mode
        run_complete_analysis(
            args.bam,
            output_dir=args.output,
            keep_temp_files=args.keep_temp,
            reuse_existing=args.reuse,
            allow_multiple=args.multiple
        )
    else:
        # Batch mode
        run_batch_analysis(
            args.dir,
            output_parent_dir=args.output,
            pattern=args.pattern,
            keep_temp_files=args.keep_temp,
            reuse_existing=args.reuse,
            allow_multiple=args.multiple
        )

if __name__ == "__main__":
    main()
