from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as build_ext_orig
from setuptools.command.install import install
import subprocess
import os
import glob
import shutil
import re

def get_version():
    version_file = os.path.join(os.path.dirname(__file__), 'aindex', '__init__.py')
    with open(version_file, 'r') as f:
        version_content = f.read()
    # Используем регулярное выражение для поиска строки с версией
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_content, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

class build_ext(build_ext_orig):
    def run(self):
        subprocess.check_call(['make', 'clean'])
        subprocess.check_call(['make', 'all'])
        build_lib = self.build_lib
        package_dir = os.path.join(build_lib, 'aindex', 'core')
        os.makedirs(package_dir, exist_ok=True)
        so_files = glob.glob(os.path.join('aindex', 'core', 'python_wrapper*.so'))
        if so_files:
            shutil.copy(so_files[0], os.path.join(package_dir, 'python_wrapper.so'))

class CustomInstall(install):
    def run(self):
        install.run(self)
        bin_dir = os.path.join(self.install_scripts, 'bin')
        os.makedirs(bin_dir, exist_ok=True)
        for file in glob.glob('bin/*'):
            shutil.copy(file, bin_dir)

setup(
    name='aindex2',
    version=get_version(),
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
        '': ['bin/*'],
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
        ],
    },
)