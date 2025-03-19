# r_repeat/unit_analysis package

# Import all functions from unit_analysis module to make them available at package level
from .unit_analysis import (
    create_motif_reference, 
    align_to_motif,
    process_motif_alignments, 
    correspond_units_to_regions,
    generate_unit_count_visualization,
    generate_unit_length_visualization,
    generate_region_length_visualization,
    find_long_units
)

__all__ = [
    'create_motif_reference',
    'align_to_motif',
    'process_motif_alignments',
    'correspond_units_to_regions',
    'generate_unit_count_visualization',
    'generate_unit_length_visualization',
    'generate_region_length_visualization',
    'find_long_units'
]
