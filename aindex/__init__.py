# aindex/__init__.py

# Import the main AIndex class
try:
    from .core.aindex import AIndex
except ImportError as e:
    # Re-raise the import error with helpful message
    raise ImportError(
        f"Failed to import AIndex: {e}\n"
        "This package requires Linux or macOS. "
        "Windows support is planned for future releases. "
        "For Windows users, please use WSL or Docker."
    ) from e

__version__ = '1.4.3'
__all__ = ['AIndex']