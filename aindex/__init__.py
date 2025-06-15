# aindex/__init__.py
import platform
import os

# Try to import the main AIndex class with Windows fallback
try:
    from .core.aindex import AIndex
except ImportError as e:
    # Check if we're on Windows (which has Python-only build)
    is_windows = platform.system() == 'Windows'
    
    if is_windows:
        print("Windows detected: C++ extensions not available. Using Python-only functionality.")
        print("For full functionality, please use Linux or macOS.")
        from .core.aindex_cpp_fallback import AIndex
    else:
        # Re-raise the original import error for other cases
        raise e

__version__ = '1.4.0'
__all__ = ['AIndex']