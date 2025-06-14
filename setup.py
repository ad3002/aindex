#!/usr/bin/env python
"""Setup script for aindex2 package."""

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as build_ext_orig
from setuptools.command.install import install
import subprocess
import os
import glob
import shutil
import sys

def check_dependencies():
    """Check if required build dependencies are available"""
    missing_deps = []
    
    for cmd in ['make', 'cmake', 'g++', 'git']:
        try:
            subprocess.check_call([cmd, '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing_deps.append(cmd)
    
    return missing_deps

def install_colab_dependencies():
    """Install missing dependencies in Google Colab environment"""
    print("Detected Google Colab environment. Installing build dependencies...")
    
    try:
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
        
        # Check if we're in a cibuildwheel environment
        in_cibw = os.environ.get('CIBUILDWHEEL', '0') == '1'
        
        try:
            if in_cibw:
                print("Building in cibuildwheel environment - building only pybind11 extension")
                subprocess.check_call(['make', 'clean'])
                subprocess.check_call(['make', 'pybind11'])
            else:
                subprocess.check_call(['make', 'clean'])
                subprocess.check_call(['make', 'all'])
        except subprocess.CalledProcessError as e:
            print(f"Build failed with error: {e}")
            print("Attempting to build with verbose output...")
            try:
                subprocess.check_call(['make', 'clean'])
                if in_cibw:
                    subprocess.check_call(['make', 'pybind11', 'VERBOSE=1'])
                else:
                    subprocess.check_call(['make', 'all', 'VERBOSE=1'])
            except subprocess.CalledProcessError as e2:
                raise RuntimeError(f"Failed to build C++ extensions: {e2}")
        
        # CRITICAL: Call parent's run() method for metadata generation
        super().run()
        
        # Copy compiled files
        self._copy_built_files()
    
    def _copy_built_files(self):
        """Copy built files to the appropriate locations"""
        build_lib = self.build_lib
        package_dir = os.path.join(build_lib, 'aindex', 'core')
        os.makedirs(package_dir, exist_ok=True)
        
        # Copy the pybind11 extension
        pybind11_files = glob.glob(os.path.join('aindex', 'core', 'aindex_cpp*.so'))
        if pybind11_files:
            for f in pybind11_files:
                shutil.copy(f, os.path.join(package_dir, os.path.basename(f)))
                print(f"Copied pybind11 extension: {f}")
        
        # Copy binaries to package
        pkg_bin_dir = os.path.join(build_lib, 'aindex', 'bin')
        os.makedirs(pkg_bin_dir, exist_ok=True)
        
        if os.path.exists('bin'):
            for file in glob.glob('bin/*'):
                dest_file = os.path.join(pkg_bin_dir, os.path.basename(file))
                shutil.copy2(file, dest_file)
                print(f"Copied binary: {os.path.basename(file)}")

class CustomInstall(install):
    def run(self):
        install.run(self)
        # Copy bin files to package data directory
        pkg_bin_dir = os.path.join(self.install_lib, 'aindex', 'bin')
        os.makedirs(pkg_bin_dir, exist_ok=True)
        
        if os.path.exists('bin'):
            for file in glob.glob('bin/*'):
                dest_file = os.path.join(pkg_bin_dir, os.path.basename(file))
                shutil.copy2(file, dest_file)
                os.chmod(dest_file, 0o755)
                print(f"Installed binary: {dest_file}")

# Minimal setup() call - most configuration comes from setup.cfg
setup(
    ext_modules=[
        Extension(
            'aindex.core.aindex_cpp', 
            sources=[],  # Built by Makefile
        ),
    ],
    cmdclass={
        'build_ext': build_ext,
        'install': CustomInstall,
    },
)