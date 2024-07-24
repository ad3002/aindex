from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
import subprocess
import os

class build_ext(build_ext_orig):
    def build_extensions(self):
        # Run makefile commands to build C++ components
        subprocess.check_call(['make', '-C', 'external/emphf'])
        subprocess.check_call(['make'])
        super().build_extensions()

# Define your C++ extensions
extensions = [
    Extension(
        'aindex.core.python_wrapper',
        sources=['src/python_wrapper.cpp'],
        include_dirs=['src', 'external/emphf'],
        language='c++', 
        extra_compile_args=['-std=c++17', '-O3', '-Wall'],
    ),
]

setup(
    name='aindex2',
    version='1.0.1',
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
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'Programming Language :: Python',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
