"""
Batch r-repeat analysis script for processing multiple samples from a directory structure.

This script:
1. Locates samples in a directory structure like /somedirs/data/{individual_sample_numbers}/
2. Concatenates multiple BAM files belonging to the same sample
3. Runs the r-repeat analysis workflow on each concatenated sample
4. Outputs results to a directory structure where each sample has its own directory

Usage: python batch_r_repeat_analyze.py --input /somedirs/data/ --output /path/to/output
"""

import os
import glob
import argparse
import subprocess
import logging
import shutil
import tempfile
from pathlib import Path
import concurrent.futures
import sys

# Import the r_repeat_workflow module
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import r_repeat_workflow

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler("logs/batch_concat_analyze.log", mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def find_samples(input_dir):
    """
    Find all sample directories in the given input directory.
    
    Parameters:
        input_dir: Root directory containing sample subdirectories
        
    Returns:
        dict: Dictionary mapping sample names to lists of BAM files
    """
    input_dir = Path(input_dir)
    samples = {}
    
    # Search for sample directories (immediate subdirectories of input_dir)
    for sample_dir in input_dir.glob("*/"):
        if not sample_dir.is_dir():
            continue
        
        sample_name = sample_dir.name
        
        # Find all BAM files in this sample directory
        bam_files = list(sample_dir.glob("**/*.bam"))
        
        if bam_files:
            samples[sample_name] = sorted(bam_files)
            logger.info(f"Found sample {sample_name} with {len(bam_files)} BAM files")
    
    logger.info(f"Found {len(samples)} samples in total")
    return samples

def concatenate_bam_files(bam_files, output_bam, temp_dir=None, threads=4):
    """
    Concatenate multiple BAM files into a single BAM file.
    
    Parameters:
        bam_files: List of BAM files to concatenate
        output_bam: Path to output concatenated BAM file
        temp_dir: Directory for temporary files
        threads: Number of threads to use for sorting
        
    Returns:
        str: Path to the concatenated BAM file
    """
    output_bam = Path(output_bam)
    
    if len(bam_files) == 1:
        # Only one BAM file, just copy or symlink it
        if output_bam.exists():
            logger.info(f"Using existing BAM file: {output_bam}")
            return str(output_bam)
        
        logger.info(f"Only one BAM file for this sample, creating copy at: {output_bam}")
        shutil.copy2(str(bam_files[0]), str(output_bam))
        return str(output_bam)
    
    # Check if output already exists and has content
    if output_bam.exists() and output_bam.stat().st_size > 0:
        logger.info(f"Using existing concatenated BAM file: {output_bam}")
        return str(output_bam)
    
    logger.info(f"Concatenating {len(bam_files)} BAM files into {output_bam}")
    
    # Create a temporary file for the merged BAM
    with tempfile.NamedTemporaryFile(suffix=".bam", dir=temp_dir, delete=False) as temp_file:
        temp_bam = temp_file.name
    
    try:
        # Use samtools to merge BAM files
        cmd = ["samtools", "merge", "-@", str(threads), temp_bam] + [str(bam) for bam in bam_files]
        logger.info(f"Running command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Sort the merged BAM file
        sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(output_bam), temp_bam]
        logger.info(f"Running command: {' '.join(sort_cmd)}")
        
        sort_result = subprocess.run(sort_cmd, check=True, capture_output=True, text=True)
        
        # Index the sorted BAM file
        index_cmd = ["samtools", "index", str(output_bam)]
        subprocess.run(index_cmd, check=True, capture_output=True, text=True)
        
        logger.info(f"Successfully created concatenated BAM file: {output_bam}")
        return str(output_bam)
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error concatenating BAM files: {e}")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        raise
    finally:
        # Clean up the temporary BAM file
        if os.path.exists(temp_bam):
            os.remove(temp_bam)

def process_sample(sample_name, bam_files, output_dir, temp_dir, threads=4, keep_temp=False):
    """
    Process a single sample by concatenating BAM files and running r-repeat analysis.
    
    Parameters:
        sample_name: Name of the sample
        bam_files: List of BAM files for this sample
        output_dir: Directory for output files
        temp_dir: Directory for temporary files
        threads: Number of threads to use
        keep_temp: Whether to keep temporary files
        
    Returns:
        dict: Results from r-repeat analysis
    """
    try:
        # Create sample-specific directories - ensure each sample has its own directory
        sample_output_dir = Path(output_dir) / sample_name
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        sample_temp_dir = Path(temp_dir) / sample_name
        sample_temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Create path for concatenated BAM file
        concat_bam = sample_temp_dir / f"{sample_name}_concatenated.bam"
        
        # Concatenate BAM files
        concatenated_bam = concatenate_bam_files(
            bam_files, 
            concat_bam,
            temp_dir=sample_temp_dir,
            threads=threads
        )
        
        # Run r-repeat analysis - explicitly ensure output is in the sample's directory
        logger.info(f"Starting r-repeat analysis for sample {sample_name}")
        results = r_repeat_workflow.run_complete_analysis(
            concatenated_bam,
            output_dir=sample_output_dir,  # Explicitly set output to sample's directory
            keep_temp_files=keep_temp,
            reuse_existing=True,  # Let's reuse existing files for efficiency
            allow_multiple=True   # Allow multiple r-repeat regions per read
        )
        
        # Add sample name to results for easier identification
        results["sample_name"] = sample_name
        
        logger.info(f"Completed r-repeat analysis for sample {sample_name}. Results in {sample_output_dir}")
        
        # Clean up concatenated BAM if requested and it's not the only BAM file
        if not keep_temp and len(bam_files) > 1:
            logger.info(f"Removing concatenated BAM file: {concat_bam}")
            try:
                os.remove(concat_bam)
                index_file = str(concat_bam) + ".bai"
                if os.path.exists(index_file):
                    os.remove(index_file)
            except Exception as e:
                logger.warning(f"Error removing temporary BAM file: {e}")
        
        return results
    
    except Exception as e:
        logger.error(f"Error processing sample {sample_name}: {e}", exc_info=True)
        return {"status": "failed", "error": str(e), "sample_name": sample_name}

def batch_process_samples(input_dir, output_dir, temp_dir=None, threads=4, max_workers=None, 
                         keep_temp=False, parallel=True):
    """
    Batch process all samples found in the input directory.
    
    Parameters:
        input_dir: Directory containing sample subdirectories
        output_dir: Directory for output files - each sample will have its own subdirectory
        temp_dir: Directory for temporary files
        threads: Number of threads per sample processing
        max_workers: Maximum number of concurrent samples to process
        keep_temp: Whether to keep temporary files
        parallel: Whether to process samples in parallel
        
    Returns:
        dict: Results for each sample
    """
    # Create parent output and temporary directories
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a reports directory to hold cross-sample analysis
    reports_dir = output_dir / "r_repeat_reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    
    if temp_dir is None:
        temp_dir = output_dir / "temp"
    
    temp_dir = Path(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    # Create logs directory if it doesn't exist
    Path("logs").mkdir(exist_ok=True)
    
    # Find all samples
    samples = find_samples(input_dir)
    
    if not samples:
        logger.error(f"No samples found in {input_dir}")
        return {"status": "error", "message": "No samples found"}
    
    # Process all samples
    results = {}
    
    # Determine whether to process in parallel or sequentially
    if parallel and max_workers is not None and max_workers > 1:
        logger.info(f"Processing {len(samples)} samples in parallel with up to {max_workers} workers")
        
        # Process samples in parallel using ThreadPoolExecutor
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_sample = {
                executor.submit(
                    process_sample, 
                    sample_name, 
                    bam_files, 
                    output_dir, 
                    temp_dir, 
                    threads,
                    keep_temp
                ): sample_name for sample_name, bam_files in samples.items()
            }
            
            for future in concurrent.futures.as_completed(future_to_sample):
                sample_name = future_to_sample[future]
                try:
                    results[sample_name] = future.result()
                except Exception as e:
                    logger.error(f"Error processing sample {sample_name}: {e}", exc_info=True)
                    results[sample_name] = {"status": "failed", "error": str(e)}
    else:
        logger.info(f"Processing {len(samples)} samples sequentially")
        
        # Process samples sequentially
        for sample_name, bam_files in samples.items():
            results[sample_name] = process_sample(
                sample_name, 
                bam_files, 
                output_dir, 
                temp_dir, 
                threads,
                keep_temp
            )
    
    # Create a cross-sample analysis report
    cross_sample_results = {
        "total_samples": len(samples),
        "successful_samples": sum(1 for result in results.values() if result.get("status") == "success"),
        "failed_samples": sum(1 for result in results.values() if result.get("status") != "success"),
        "sample_results": {
            name: {
                "status": result.get("status", "unknown"),
                "output_directory": str(output_dir / name)  # Include explicit path to sample directory
            } 
            for name, result in results.items()
        }
    }
    
    # Write the cross-sample report to the reports directory
    report_file = reports_dir / "batch_r_repeat_analysis_report.json"
    import json
    with open(report_file, 'w') as f:
        json.dump(cross_sample_results, f, indent=2)
    
    # Create symlink to the output directory for convenience
    latest_link = Path("batch_results_latest")
    if latest_link.is_symlink() or latest_link.exists():
        latest_link.unlink()
    
    try:
        latest_link.symlink_to(output_dir.resolve())
    except:
        logger.warning("Could not create 'batch_results_latest' symlink. This may be due to filesystem limitations.")
    
    # Run batch cross-sample analysis for upstream and downstream regions
    # and save results in the reports directory
    try:
        logger.info("Starting batch cross-sample analysis for upstream regions")
        cross_sample_analysis_script = Path(__file__).parent / "r_repeat_batch_analyze_upstream.py"
        if cross_sample_analysis_script.exists():
            cmd = ["python", str(cross_sample_analysis_script), "--input", str(output_dir), "--output", str(reports_dir / "upstream")]
            subprocess.run(cmd, check=True)
    except Exception as e:
        logger.error(f"Error running batch upstream analysis: {e}")
    
    try:
        logger.info("Starting batch cross-sample analysis for downstream regions")
        cross_sample_analysis_script = Path(__file__).parent / "r_repeat_batch_analyze_downstream.py"
        if cross_sample_analysis_script.exists():
            cmd = ["python", str(cross_sample_analysis_script), "--input", str(output_dir), "--output", str(reports_dir / "downstream")]
            subprocess.run(cmd, check=True)
    except Exception as e:
        logger.error(f"Error running batch downstream analysis: {e}")
    
    logger.info(f"Batch r-repeat analysis complete. Processed {len(samples)} samples.")
    logger.info(f"Successful samples: {cross_sample_results['successful_samples']}")
    logger.info(f"Failed samples: {cross_sample_results['failed_samples']}")
    logger.info(f"Results for each sample are in separate directories under {output_dir}")
    logger.info(f"Cross-sample reports saved to {reports_dir}")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description="Batch process samples by concatenating BAMs and running r-repeat analysis"
    )
    
    parser.add_argument("--input", "-i", required=True,
                      help="Directory containing sample subdirectories with BAM files")
    
    # Add description for the output directory to clarify each sample gets its own directory
    parser.add_argument("--output", "-o", required=True,
                      help="Parent output directory - each sample will get a subdirectory with its name")
    
    parser.add_argument("--temp", "-t",
                      help="Directory for temporary files (default: output_dir/temp)")
    
    parser.add_argument("--threads", type=int, default=4,
                      help="Number of threads to use per sample (default: 4)")
    
    parser.add_argument("--max-workers", type=int, default=1,
                      help="Maximum number of samples to process in parallel (default: 1)")
    
    parser.add_argument("--keep-temp", action="store_true",
                      help="Keep temporary files after analysis")
    
    parser.add_argument("--sequential", action="store_true",
                      help="Process samples sequentially even if max-workers > 1")
    
    args = parser.parse_args()
    
    batch_process_samples(
        args.input,
        args.output,
        args.temp,
        args.threads,
        args.max_workers,
        args.keep_temp,
        not args.sequential  # Process in parallel if not sequential
    )

if __name__ == "__main__":
    main()
