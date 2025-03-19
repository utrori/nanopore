# R-repeat analysis package
"""
R-repeat analysis package for Oxford Nanopore sequencing data.
This package provides tools for detecting, analyzing, and visualizing r-repeat regions in nanopore reads.
"""

# Import core modules for easier access
from .core.workflow import run_complete_analysis, run_batch_analysis
from .detection.detection import detect_r_repeats
from .methylation.methylation_analysis import analyze_methylation

# Package information
__version__ = '1.0.0'
__author__ = 'Your Name'
__email__ = 'your.email@example.com'
