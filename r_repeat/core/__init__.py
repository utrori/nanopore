# r_repeat/core package

# To avoid circular imports, we'll define what should be imported 
# but won't actually import them here

__all__ = ['run_complete_analysis', 'run_batch_analysis']

# Use __getattr__ to lazy load these attributes
def __getattr__(name):
    if name in __all__:
        from .workflow import run_complete_analysis, run_batch_analysis
        if name == 'run_complete_analysis':
            return run_complete_analysis
        elif name == 'run_batch_analysis':
            return run_batch_analysis
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
