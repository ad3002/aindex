# Google Colab Installation Guide

To install aindex in Google Colab, follow these steps:

## ⚠️ IMPORTANT: cmake Conflict in Google Colab

Google Colab has a conflict between the Python `cmake` package and system cmake. This needs to be fixed first:

```python
# Quick fix for cmake conflict
!wget https://raw.githubusercontent.com/ad3002/aindex/main/fix_cmake_colab.py
!python fix_cmake_colab.py
```

## Method 1: Automatic Installation (Recommended)

```python
# Download and run installation script
!wget https://raw.githubusercontent.com/ad3002/aindex/main/install_colab.py
!python install_colab.py
```

## Method 2: Manual Installation

```python
# 1. Fix cmake conflict
!pip uninstall -y cmake || true
!apt-get update
!apt-get install -y build-essential cmake git python3-dev

# 2. Clone repository
!git clone https://github.com/ad3002/aindex.git
%cd aindex

# 3. Install Python dependencies
!pip install -r requirements.txt

# 4. Build package
!make clean
!make all

# 5. Install package
!pip install -e .
```

## Method 3: Installation from Source Code

If you have a local copy:

```python
# Upload files to Google Colab (use Files -> Upload)
# Then:
!pip uninstall -y cmake || true
!apt-get update
!apt-get install -y build-essential cmake git python3-dev
!pip install -r requirements.txt
!make clean && make all
!pip install -e .
```

## Installation Verification

```python
import aindex
print("aindex successfully installed!")
```

## Troubleshooting

### 1. Error "Failed to clone emphf repository"
- Check your internet connection
- Try running: `!git config --global http.timeout 300`

### 2. C++ Build Error
- Make sure build-essential and cmake are installed
- Try: `!apt-get install -y g++ cmake make`

### 3. Python Wrapper Error
- Verify that the build completed successfully
- Make sure the file `aindex/core/python_wrapper.so` is created

### 4. Import Error
- Restart runtime: Runtime -> Restart Runtime
- Reinstall package: `!pip install -e . --force-reinstall`
