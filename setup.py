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

def get_version():
    version_file = os.path.join(os.path.dirname(__file__), 'aindex', '__init__.py')
    with open(version_file, 'r') as f:
        version_content = f.read()
    # Используем регулярное выражение для поиска строки с версией
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_content, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

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
        
        try:
            subprocess.check_call(['make', 'clean'])
            subprocess.check_call(['make', 'all'])  # Build everything including binaries
        except subprocess.CalledProcessError as e:
            print(f"Build failed with error: {e}")
            print("Attempting to build with verbose output...")
            try:
                subprocess.check_call(['make', 'clean'])
                subprocess.check_call(['make', 'all', 'VERBOSE=1'])
            except subprocess.CalledProcessError as e2:
                raise RuntimeError(f"Failed to build C++ extensions: {e2}")
        
        build_lib = self.build_lib
        package_dir = os.path.join(build_lib, 'aindex', 'core')
        os.makedirs(package_dir, exist_ok=True)
        
        # Copy the pybind11 extension (modern API)
        pybind11_files = glob.glob(os.path.join('aindex', 'core', 'aindex_cpp*.so'))
        if pybind11_files:
            shutil.copy(pybind11_files[0], os.path.join(package_dir, os.path.basename(pybind11_files[0])))
            print(f"Copied pybind11 extension: {pybind11_files[0]}")
        
        # Copy binaries to package
        pkg_bin_dir = os.path.join(build_lib, 'aindex', 'bin')
        os.makedirs(pkg_bin_dir, exist_ok=True)
        
        if os.path.exists('bin'):
            for file in glob.glob('bin/*'):
                dest_file = os.path.join(pkg_bin_dir, os.path.basename(file))
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

setup(
    ext_modules=[
        Extension('aindex.core.aindex_cpp', sources=[]),  # Built by Makefile
    ],
    cmdclass={
        'build_ext': build_ext,
        'install': CustomInstall,
    },
    include_package_data=True,
    package_data={
        'aindex.core': ['*.so', 'aindex_cpp*.so'],
        'aindex': ['bin/*'],
    },
)