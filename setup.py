from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as build_ext_orig
import subprocess
import os
import shutil
import glob
from setuptools.command.install import install


class build_ext(build_ext_orig):
    def run(self):
        # Run makefile commands to build C++ components
        subprocess.check_call(['make', 'clean'])
        subprocess.check_call(['make', 'all'])
        
        # Copy the pre-built shared object to the correct location
        build_lib = self.build_lib
        package_dir = os.path.join(build_lib, 'aindex', 'core')
        os.makedirs(package_dir, exist_ok=True)
        
        # Find the .so file and rename it
        so_file = glob.glob(os.path.join('aindex', 'core', 'python_wrapper*.so'))[0]
        shutil.copy(so_file, os.path.join(package_dir, 'python_wrapper.so'))

class CustomInstall(install):
    def run(self):
        install.run(self)
        # Ensure bin directory exists in the installation
        bin_dir = os.path.join(self.install_scripts, 'bin')
        os.makedirs(bin_dir, exist_ok=True)
        
        # Copy bin files
        for file in glob.glob('bin/*'):
            shutil.copy(file, bin_dir)

setup(
    name='aindex2',
    version='1.0.4',
    description='Perfect hash based index for genome data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Aleksey Komissarov',
    author_email='ad3002@gmail.com',
    url='https://github.com/ad3002/aindex',
    packages=find_packages(),
    ext_modules=[Extension('aindex.core.python_wrapper', sources=[])],
    cmdclass={
        'build_ext': build_ext,
        'install': CustomInstall,
    },
    install_requires=open('requirements.txt').read().splitlines(),
    include_package_data=True,
    package_data={
        'aindex.core': ['*.so'],
    },
    data_files=[
        ('bin', glob.glob('bin/*')),
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
    ],
    entry_points={
        'console_scripts': [
            'compute_index_py=scripts.compute_index:main',
            'compute_aindex_py=scripts.compute_aindex:main',
            'reads_to_fasta=scripts.reads_to_fasta:main',
        ],
    },
)