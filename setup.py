from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
import subprocess
import os
import shutil

class build_ext(build_ext_orig):
    def run(self):
        # Run makefile commands to build C++ components
        subprocess.check_call(['make', 'clean'])
        subprocess.check_call(['make', 'all'])
        
        # Copy the pre-built shared object to the correct location
        build_lib = self.build_lib
        package_dir = os.path.join(build_lib, 'aindex', 'core')
        os.makedirs(package_dir, exist_ok=True)
        shutil.copy('aindex/core/python_wrapper.so', package_dir)

# Define your C++ extensions, but don't actually build them
extensions = [
    Extension(
        'aindex.core.python_wrapper',
        sources=[],  # Empty sources as we're not building it here
    ),
]

setup(
    name='aindex2',
    version='1.0.3',
    description='Perfect hash based index for genome data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Aleksey Komissarov',
    author_email='ad3002@gmail.com',
    url='https://github.com/ad3002/aindex',
    packages=['aindex', 'aindex.core'],
    ext_modules=extensions,
    cmdclass={'build_ext': build_ext},
    install_requires=open('requirements.txt').read().splitlines(),
    include_package_data=True,
    package_data={
        'aindex.core': ['*.so'],
    },
    data_files=[
        ('bin', [
            'bin/compute_index.exe', 
            'bin/compute_aindex.exe', 
            'bin/compute_reads.exe',
            'bin/compute_aindex.py', 
            'bin/compute_index.py', 
            'bin/reads_to_fasta.py'
        ]),
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
    ],
    entry_points={
        'console_scripts': [
            'compute_index=scripts.compute_index:main',
            'compute_aindex=scripts.compute_aindex:main',
            'reads_to_fasta=scripts.reads_to_fasta:main',
        ],
    },
)