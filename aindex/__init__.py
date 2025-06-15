# aindex/__init__.py
import platform
import os

# Try to import the main AIndex class with Windows fallback
try:
    from .core.aindex import AIndex
except ImportError as e:
    # Check if we're on Windows with Python-only build
    is_windows = platform.system() == 'Windows'
    windows_python_only = os.environ.get('WINDOWS_PYTHON_ONLY', '0') == '1'
    
    if is_windows and windows_python_only:
        print("Windows Python-only build detected - C++ extensions not available")
        from .core.aindex_cpp_fallback import AIndex
    else:
        # Re-raise the original import error for other cases
        raise e

__version__ = '1.4.0'
__all__ = ['AIndex']