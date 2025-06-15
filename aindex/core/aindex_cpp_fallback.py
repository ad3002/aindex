# Windows fallback module for aindex_cpp
# This module provides minimal functionality when C++ extensions are not available

class WindowsFallback:
    """Fallback implementation for Windows when C++ extensions are not built"""
    
    def __init__(self):
        raise NotImplementedError(
            "C++ extensions are not available on Windows. "
            "This Windows build provides Python-only functionality. "
            "For full functionality, please use Linux or macOS builds."
        )

# For compatibility, export the same interface
AIndex = WindowsFallback
