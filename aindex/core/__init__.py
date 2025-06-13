# aindex/core/__init__.py
# Import only the main classes and functions, not everything
from .aindex import AIndex, get_revcomp, hamming_distance, Strand

__all__ = ['AIndex', 'get_revcomp', 'hamming_distance', 'Strand']