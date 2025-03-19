"""
R-repeat analysis package for Oxford Nanopore sequencing data.
This package provides tools for detecting, analyzing, and visualizing r-repeat regions in nanopore reads.
"""

# To avoid circular imports, we'll define what should be imported directly
# But won't actually import them here - they'll be imported when needed

__all__ = [
    'run_complete_analysis',
    'run_batch_analysis', 
    'detect_r_repeats',
    'analyze_methylation'
]

# Package information
__version__ = '1.0.0'
__author__ = 'Your Name'
__email__ = 'your.email@example.com'

# Define a function to load the package components when needed
def _load_components():
    from r_repeat.core.workflow import run_complete_analysis, run_batch_analysis
    from r_repeat.detection.detection import detect_r_repeats
    from r_repeat.methylation.methylation_analysis import analyze_methylation
    
    return {
        'run_complete_analysis': run_complete_analysis,
        'run_batch_analysis': run_batch_analysis,
        'detect_r_repeats': detect_r_repeats,
        'analyze_methylation': analyze_methylation
    }

# Load components only when directly accessed
def __getattr__(name):
    if name in __all__:
        components = _load_components()
        if name in components:
            return components[name]
    raise AttributeError(f"module {__name__} has no attribute {name}")
