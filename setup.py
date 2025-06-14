from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as build_ext_orig
from setuptools.command.install import install
import subprocess
import os
import glob
import shutil
import re
import sys
import platform

# Important: define package metadata at the beginning of the file
PACKAGE_NAME = "aindex2"
PACKAGE_VERSION = "1.3.22"

def check_dependencies():
    """Check if required build dependencies are available"""
    missing_deps = []
    
    # Check for make
    try:
        subprocess.check_call(['make', '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing_deps.append('make')
    
    # Check for cmake
    try:
        subprocess.check_call(['cmake', '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing_deps.append('cmake')
    
    # Check for g++
    try:
        subprocess.check_call(['g++', '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing_deps.append('g++')
    
    # Check for git
    try:
        subprocess.check_call(['git', '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing_deps.append('git')
    
    return missing_deps

def install_colab_dependencies():
    """Install missing dependencies in Google Colab environment"""
    print("Detected Google Colab environment. Installing build dependencies...")
    
    try:
        # Install build essentials
        subprocess.check_call(['apt-get', 'update'], stdout=subprocess.DEVNULL)
        subprocess.check_call(['apt-get', 'install', '-y', 'build-essential', 'cmake', 'git'], 
                            stdout=subprocess.DEVNULL)
        print("Build dependencies installed successfully.")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Failed to install dependencies: {e}")
        return False

class build_ext(build_ext_orig):
    def run(self):
        # Check if we're in Google Colab
        in_colab = 'google.colab' in sys.modules
        
        if in_colab:
            print("Google Colab environment detected.")
            missing_deps = check_dependencies()
            if missing_deps:
                print(f"Missing dependencies: {', '.join(missing_deps)}")
                if not install_colab_dependencies():
                    raise RuntimeError("Failed to install required build dependencies")
        
        # Check if we're in a cibuildwheel environment or if binaries already exist
        in_cibw = os.environ.get('CIBUILDWHEEL', '0') == '1'
        
        try:
            if in_cibw:
                # In cibuildwheel, external dependencies should be built in BEFORE_ALL
                # We only need to build the pybind11 extension
                print("Building in cibuildwheel environment - building only pybind11 extension")
                subprocess.check_call(['make', 'clean'])
                subprocess.check_call(['make', 'pybind11'])
            else:
                # Regular build - build everything
                subprocess.check_call(['make', 'clean'])
                subprocess.check_call(['make'])  # Simplified make call
        except subprocess.CalledProcessError as e:
            print(f"Build failed with error: {e}")
            print("Attempting to build with verbose output...")
            try:
                subprocess.check_call(['make', 'clean'])
                if in_cibw:
                    subprocess.check_call(['make', 'pybind11', 'VERBOSE=1'])
                else:
                    subprocess.check_call(['make', 'VERBOSE=1'])  # Simplified make call
            except subprocess.CalledProcessError as e2:
                raise RuntimeError(f"Failed to build C++ extensions: {e2}")
        
        # Important: call the parent run() method for correct metadata handling
        super().run()
        
        build_lib = self.build_lib
        package_dir = os.path.join(build_lib, 'aindex', 'core')
        os.makedirs(package_dir, exist_ok=True)
        
        # Copy the pybind11 extension (modern API)
        pybind11_files = glob.glob(os.path.join('aindex', 'core', 'aindex_cpp*.so'))
        if pybind11_files:
            shutil.copy(pybind11_files[0], os.path.join(package_dir, os.path.basename(pybind11_files[0])))
            print(f"Copied pybind11 extension: {pybind11_files[0]}")
        
        # Copy binaries to package bin directory (important!)
        pkg_bin_dir = os.path.join(build_lib, 'aindex', 'bin')
        os.makedirs(pkg_bin_dir, exist_ok=True)
        
        # Try multiple sources for binaries
        binary_sources = ['bin', 'aindex/bin']
        binaries_copied = 0
        
        for source_dir in binary_sources:
            if os.path.exists(source_dir):
                for file in glob.glob(os.path.join(source_dir, '*')):
                    if os.path.isfile(file):  # Only copy files, not directories
                        dest_file = os.path.join(pkg_bin_dir, os.path.basename(file))
                        shutil.copy2(file, dest_file)
                        # Make executable on Unix-like systems
                        if not file.endswith('.py'):
                            os.chmod(dest_file, 0o755)
                        print(f"Copied binary to package: {os.path.basename(file)} from {source_dir}")
                        binaries_copied += 1
        
        if binaries_copied == 0:
            print("Warning: No binaries found to copy from bin/ or aindex/bin/")
            # In CI, try to ensure binaries exist in aindex/bin
            if in_cibw or os.environ.get('CI'):
                print("CI environment detected, ensuring binaries are copied to aindex/bin")
                if os.path.exists('bin'):
                    # Copy from bin to aindex/bin and then to package
                    os.makedirs('aindex/bin', exist_ok=True)
                    for file in glob.glob('bin/*'):
                        if os.path.isfile(file):
                            dest_in_source = os.path.join('aindex/bin', os.path.basename(file))
                            shutil.copy2(file, dest_in_source)
                            # Now copy to build directory
                            dest_file = os.path.join(pkg_bin_dir, os.path.basename(file))
                            shutil.copy2(file, dest_file)
                            if not file.endswith('.py'):
                                os.chmod(dest_file, 0o755)
                            print(f"CI: Copied binary {os.path.basename(file)}")
                            binaries_copied += 1
        
        # Old method - copy to separate pkg_bin_dir (kept for compatibility)
        pkg_bin_dir_old = os.path.join(build_lib, 'aindex', 'bin')
        os.makedirs(pkg_bin_dir_old, exist_ok=True)
        
        if os.path.exists('bin'):
            for file in glob.glob('bin/*'):
                dest_file = os.path.join(pkg_bin_dir_old, os.path.basename(file))
                shutil.copy2(file, dest_file)
                print(f"Copied binary: {os.path.basename(file)}")
        else:
            print("Warning: No pybind11 extension found.")

class CustomInstall(install):
    def run(self):
        install.run(self)
        # Copy bin files to package data directory
        pkg_bin_dir = os.path.join(self.install_lib, 'aindex', 'bin')
        os.makedirs(pkg_bin_dir, exist_ok=True)
        
        # Copy all binaries
        if os.path.exists('bin'):
            for file in glob.glob('bin/*'):
                dest_file = os.path.join(pkg_bin_dir, os.path.basename(file))
                shutil.copy2(file, dest_file)
                # Make executable
                os.chmod(dest_file, 0o755)
                print(f"Installed binary: {dest_file}")

# Read README for long_description
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name=PACKAGE_NAME,
    version=PACKAGE_VERSION,
    description="Perfect hash based index for genome data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Aleksey Komissarov",
    author_email="ad3002@gmail.com",
    url="https://github.com/ad3002/aindex",
    packages=find_packages(),
    ext_modules=[
        Extension(
            'aindex.core.aindex_cpp', 
            sources=[],  # Built by Makefile
            # Add empty include_dirs and library_dirs for compatibility
            include_dirs=[],
            library_dirs=[],
        ),
    ],
    cmdclass={
        'build_ext': build_ext,
        'install': CustomInstall,
    },
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "intervaltree==3.1.0", 
        "editdistance==0.8.1",
        "psutil>=5.8.0",
    ],
    include_package_data=True,
    package_data={
        'aindex.core': ['*.so', 'aindex_cpp*.so'],
        'aindex': ['bin/*', 'bin/**/*'],
    },
    entry_points={
        'console_scripts': [
            'aindex=aindex.cli:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9", 
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: C++",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
    ],
    # Add zip_safe=False for correct work with binary extensions
    zip_safe=False,
)