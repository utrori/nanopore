"""
Simple batch r-repeat analysis script for processing multiple samples.

Assumes input directory structure like:
/path/to/input/{sample_name}/{sample_bam_files...}

And creates outputs in:
/path/to/output/{sample_name}/...
"""

import os
import subprocess
import logging
import shutil
from pathlib import Path
import sys
import concurrent.futures

# Fix the import issue by making sure the parent directory is in the Python path
current_dir = Path(__file__).resolve()
project_root = current_dir.parent.parent.parent  # This gets to /home/owner/nanopore_rDNA
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Now we can import from r_repeat
try:
    from r_repeat.core.workflow import run_complete_analysis
except ModuleNotFoundError:
    # Fallback approach - try direct import using relative path if the package import fails
    sys.path.insert(0, str(current_dir.parent.parent))  # Add r_repeat directory to path
    from core.workflow import run_complete_analysis

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def process_sample(sample_dir, output_parent_dir, keep_temp=False):
    """
    Process a single sample directory containing BAM files.
    
    Args:
        sample_dir: Path to sample directory containing BAM files
        output_parent_dir: Parent directory for output
        keep_temp: Whether to keep temporary files
    """
    sample_dir = Path(sample_dir)
    sample_name = sample_dir.name
    
    # Find all BAM files for this sample
    bam_files = list(sample_dir.glob("*.bam"))
    
    if not bam_files:
        logger.warning(f"No BAM files found in {sample_dir}")
        return {"status": "error", "sample": sample_name, "message": "No BAM files found"}
    
    # If multiple BAM files, concatenate them
    if len(bam_files) > 1:
        logger.info(f"Found {len(bam_files)} BAM files for {sample_name}")
        
        # Create temporary directory for concatenated BAM
        temp_dir = Path(output_parent_dir) / sample_name / "temp"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Concatenate BAMs
        concat_bam = temp_dir / f"{sample_name}_concatenated.bam"
        
        # Use samtools merge
        cmd = ["samtools", "merge", "-o", str(concat_bam)] + [str(f) for f in bam_files]
        logger.info(f"Concatenating BAMs with command: {' '.join(cmd)}")
        
        subprocess.run(cmd, check=True)
        bam_path = concat_bam
    else:
        # Only one BAM file, use it directly
        bam_path = bam_files[0]
    
    # Set up output directory for this sample
    output_dir = Path(output_parent_dir) / sample_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run the analysis
    logger.info(f"Processing sample {sample_name}")
    results = run_complete_analysis(
        bam_path,
        output_dir=output_dir,
        keep_temp_files=keep_temp,
        reuse_existing=True,
        allow_multiple=True
    )
    
    logger.info(f"Completed analysis for sample {sample_name}")
    return {"status": "success", "sample": sample_name, "results": results}

def batch_process(input_dir, output_dir, parallel=False, max_workers=1, keep_temp=False):
    """
    Process all sample directories in the input directory.
    
    Args:
        input_dir: Directory containing sample subdirectories
        output_dir: Parent directory for output files
        parallel: Whether to process samples in parallel
        max_workers: Maximum number of concurrent processes
        keep_temp: Whether to keep temporary files
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all subdirectories (each is a sample)
    sample_dirs = [d for d in input_dir.iterdir() if d.is_dir()]
    logger.info(f"Found {len(sample_dirs)} sample directories in {input_dir}")
    
    results = {}
    
    if parallel and max_workers > 1:
        # Process samples in parallel
        logger.info(f"Processing samples in parallel with {max_workers} workers")
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_sample, sample_dir, output_dir, keep_temp): sample_dir.name
                for sample_dir in sample_dirs
            }
            
            for future in concurrent.futures.as_completed(futures):
                sample_name = futures[future]
                try:
                    results[sample_name] = future.result()
                except Exception as e:
                    logger.error(f"Error processing sample {sample_name}: {str(e)}")
                    results[sample_name] = {"status": "error", "sample": sample_name, "error": str(e)}
    else:
        # Process samples sequentially
        logger.info("Processing samples sequentially")
        for sample_dir in sample_dirs:
            sample_name = sample_dir.name
            try:
                results[sample_name] = process_sample(sample_dir, output_dir, keep_temp)
            except Exception as e:
                logger.error(f"Error processing sample {sample_name}: {str(e)}")
                results[sample_name] = {"status": "error", "sample": sample_name, "error": str(e)}
    
    logger.info(f"Completed batch processing of {len(sample_dirs)} samples")
    return results

if __name__ == "__main__":
    # Example usage with hardcoded paths
    INPUT_DIR = "/path/to/input"
    OUTPUT_DIR = "/path/to/output"
    PARALLEL = False
    MAX_WORKERS = 1
    KEEP_TEMP = False
    
    batch_process(INPUT_DIR, OUTPUT_DIR, PARALLEL, MAX_WORKERS, KEEP_TEMP)

    # Uncomment if you want command line arguments
    """
    import argparse
    parser = argparse.ArgumentParser(description="Batch process r-repeat samples")
    parser.add_argument("--input", "-i", required=True, help="Input directory containing sample directories")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--parallel", "-p", action="store_true", help="Process samples in parallel")
    parser.add_argument("--workers", "-w", type=int, default=1, help="Maximum number of workers for parallel processing")
    parser.add_argument("--keep-temp", "-k", action="store_true", help="Keep temporary files")
    
    args = parser.parse_args()
    batch_process(args.input, args.output, args.parallel, args.workers, args.keep_temp)
    """
