#!/usr/bin/env python3
"""
aindex command-line interface
Unified CLI for all aindex tools and utilities
"""

import argparse
import sys
import os
import subprocess
import shutil
from pathlib import Path
import pkg_resources
import importlib.util


def get_bin_path():
    """Get the path to the bin directory containing executables"""
    
    # Method 1: Try to find in installed package directory first (via aindex module path)
    try:
        import aindex
        package_dir = Path(aindex.__file__).parent
        bin_dir = package_dir / "bin"
        if bin_dir.exists() and any(bin_dir.iterdir()):
            return bin_dir
    except:
        pass
    
    # Method 2: Try via pkg_resources (installed package)
    try:
        package_path = pkg_resources.resource_filename('aindex', 'bin')
        bin_dir = Path(package_path)
        if bin_dir.exists() and any(bin_dir.iterdir()):
            return bin_dir
    except:
        pass
    
    # Method 3: Try to find in site-packages
    try:
        import site
        import sys
        for site_dir in site.getsitepackages() + [site.getusersitepackages()]:
            if site_dir:
                bin_dir = Path(site_dir) / "aindex" / "bin" 
                if bin_dir.exists() and any(bin_dir.iterdir()):
                    return bin_dir
    except:
        pass
    
    # Method 4: Fallback to development environment (relative to this file)
    current_dir = Path(__file__).parent.parent
    bin_dir = current_dir / "bin"
    
    if bin_dir.exists():
        return bin_dir
    
    # Fallback to current directory
    return Path.cwd() / "bin"


def run_executable(exe_name, args):
    """Run a binary executable from bin directory"""
    bin_dir = get_bin_path()
    
    # Try different executable extensions based on platform
    candidates = [exe_name]
    if exe_name.endswith('.exe'):
        # If requested with .exe, also try without extension
        candidates.append(exe_name[:-4])
    else:
        # If requested without extension, also try with .exe
        candidates.append(exe_name + '.exe')
    
    exe_path = None
    for candidate in candidates:
        candidate_path = bin_dir / candidate
        if candidate_path.exists():
            exe_path = candidate_path
            break
    
    if exe_path is None:
        print(f"Error: Executable {exe_name} not found in {bin_dir}")
        print(f"Tried: {[str(bin_dir / c) for c in candidates]}")
        return 1
    
    try:
        # Run the executable with the provided arguments
        result = subprocess.run([str(exe_path)] + args, check=False)
        return result.returncode
    except Exception as e:
        print(f"Error running {exe_path.name}: {e}")
        return 1


def run_python_script(script_name, args):
    """Run a Python script from scripts directory"""
    # Try to find scripts directory
    current_dir = Path(__file__).parent.parent
    scripts_dir = current_dir / "scripts"
    
    if not scripts_dir.exists():
        # Try in bin directory
        bin_dir = get_bin_path()
        script_path = bin_dir / script_name
    else:
        script_path = scripts_dir / script_name
    
    if not script_path.exists():
        print(f"Error: Script {script_name} not found")
        return 1
    
    try:
        # Run the Python script with the provided arguments
        result = subprocess.run([sys.executable, str(script_path)] + args, check=False)
        return result.returncode
    except Exception as e:
        print(f"Error running {script_name}: {e}")
        return 1


def cmd_compute_aindex(args):
    """Compute aindex for k-mer analysis"""
    parser = argparse.ArgumentParser(
        prog='aindex compute-aindex',
        description='Compute aindex for k-mer analysis (supports 13-mer and 23-mer modes)'
    )
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ/FASTA files (comma-separated)')
    parser.add_argument('-t', '--type', default='fastq', choices=['fastq', 'fasta'], help='Input file type')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for generated files')
    parser.add_argument('-k', '--kmer-size', type=int, choices=[13, 23], default=23, help='K-mer size (13 or 23)')
    parser.add_argument('--lu', type=int, default=2, help='Lower frequency threshold for k-mers')
    parser.add_argument('-P', '--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--use-kmer-counter', action='store_true', help='Use built-in fast k-mer counter instead of jellyfish')
    
    parsed_args = parser.parse_args(args)
    
    # Convert to arguments for the underlying script
    script_args = [
        '-i', parsed_args.input,
        '-t', parsed_args.type,
        '-o', parsed_args.output,
        '--lu', str(parsed_args.lu),
        '-P', str(parsed_args.threads)
    ]
    
    if parsed_args.use_kmer_counter:
        script_args.append('--use_kmer_counter')
    
    # Add k-mer size specific logic
    if parsed_args.kmer_size == 13:
        print(f"Computing 13-mer aindex for {parsed_args.input}")
        # For 13-mers, we might need special handling - but compute_aindex.py should handle this automatically
    else:
        print(f"Computing 23-mer aindex for {parsed_args.input}")
    
    return run_python_script('compute_aindex.py', script_args)


def cmd_compute_index(args):
    """Compute index from input data"""
    parser = argparse.ArgumentParser(
        prog='aindex compute-index',
        description='Compute index from input data'
    )
    parser.add_argument('-i', '--input', required=True, help='Input file')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    
    parsed_args = parser.parse_args(args)
    
    script_args = ['-i', parsed_args.input, '-o', parsed_args.output]
    return run_python_script('compute_index.py', script_args)


def cmd_compute_reads(args):
    """Process reads using compute_reads"""
    print("Processing reads...")
    return run_executable('compute_reads', args)


def cmd_count_kmers(args):
    """Count k-mers using fast built-in counter"""
    parser = argparse.ArgumentParser(
        prog='aindex count',
        description='Count k-mers using built-in fast counter'
    )
    parser.add_argument('-i', '--input', required=True, help='Input FASTA/FASTQ file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-k', '--kmer-size', type=int, choices=[13, 23], default=23, help='K-mer size')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads')
    
    parsed_args = parser.parse_args(args)
    
    if parsed_args.kmer_size == 13:
        print(f"Counting 13-mers in {parsed_args.input}")
        # Use specialized 13-mer counter if available
        exe_args = [parsed_args.input, parsed_args.output, str(parsed_args.threads)]
        return run_executable('count_kmers13', exe_args)
    else:
        print(f"Counting {parsed_args.kmer_size}-mers in {parsed_args.input}")
        exe_args = [parsed_args.input, parsed_args.output, str(parsed_args.kmer_size), str(parsed_args.threads)]
        return run_executable('kmer_counter', exe_args)


def cmd_build_hash(args):
    """Build perfect hash for k-mers"""
    parser = argparse.ArgumentParser(
        prog='aindex build-hash',
        description='Build perfect hash for k-mers'
    )
    parser.add_argument('-i', '--input', required=True, help='Input k-mers file')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('-k', '--kmer-size', type=int, choices=[13, 23], default=23, help='K-mer size')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    
    parsed_args = parser.parse_args(args)
    
    if parsed_args.kmer_size == 13:
        print(f"Building 13-mer hash for {parsed_args.input}")
        exe_args = [parsed_args.input, parsed_args.output, str(parsed_args.threads)]
        return run_executable('build_13mer_hash', exe_args)
    else:
        print(f"Building hash for {parsed_args.kmer_size}-mers")
        # Use general purpose hash builder
        exe_args = [parsed_args.input, parsed_args.output, str(parsed_args.threads)]
        return run_executable('compute_mphf_seq', exe_args)


def cmd_generate_kmers(args):
    """Generate all possible k-mers"""
    parser = argparse.ArgumentParser(
        prog='aindex generate',
        description='Generate all possible k-mers'
    )
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-k', '--kmer-size', type=int, choices=[13, 23], default=13, help='K-mer size')
    
    parsed_args = parser.parse_args(args)
    
    if parsed_args.kmer_size == 13:
        print(f"Generating all 13-mers to {parsed_args.output}")
        exe_args = [parsed_args.output]
        return run_executable('generate_all_13mers', exe_args)
    else:
        print(f"Generating all {parsed_args.kmer_size}-mers is not supported (too many combinations)")
        return 1


def cmd_reads_to_fasta(args):
    """Convert reads to FASTA format"""
    print("Converting reads to FASTA...")
    return run_python_script('reads_to_fasta.py', args)


def cmd_version(args):
    """Show aindex version"""
    try:
        import aindex
        print(f"aindex version {aindex.__version__}")
        
        # Show available tools
        bin_dir = get_bin_path()
        if bin_dir.exists():
            executables = [f.name for f in bin_dir.iterdir() if f.is_file() and not f.name.endswith('.py') and not f.name.startswith('__')]
            if executables:
                print(f"Available executables in {bin_dir}:")
                for exe in sorted(executables):
                    print(f"  {exe}")
        
        return 0
    except ImportError:
        print("aindex not properly installed")
        return 1


def cmd_info(args):
    """Show system and installation information"""
    try:
        import aindex
        import aindex.core.aindex_cpp as aindex_cpp
        
        print("=== aindex System Information ===")
        print(f"Version: {aindex.__version__}")
        print(f"Python: {sys.version}")
        print(f"Platform: {sys.platform}")
        
        # Test C++ module
        try:
            wrapper = aindex_cpp.AindexWrapper()
            methods = [m for m in dir(wrapper) if not m.startswith('_')]
            print(f"C++ API: Available ({len(methods)} methods)")
        except Exception as e:
            print(f"C++ API: Error - {e}")
        
        # Show detailed path information
        print("\n=== Path Information ===")
        
        # Show aindex package location
        try:
            package_dir = Path(aindex.__file__).parent
            print(f"Package location: {package_dir}")
            print(f"Package type: {'Development' if str(package_dir).endswith('workspace/aindex/aindex') else 'Installed'}")
        except:
            print("Package location: Unknown")
        
        # Show bin directory search results
        bin_dir = get_bin_path()
        print(f"Bin directory: {bin_dir}")
        print(f"Bin exists: {bin_dir.exists()}")
        
        if bin_dir.exists():
            files = list(bin_dir.iterdir())
            print(f"Bin files: {len(files)} files")
            
            # Show executables
            executables = [f.name for f in files if f.is_file() and not f.name.endswith('.py') and not f.name.startswith('__')]
            if executables:
                print("Executables:", ", ".join(sorted(executables)))
        
        return 0
    except Exception as e:
        print(f"Error getting info: {e}")
        return 1


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        prog='aindex',
        description='aindex: perfect hash based index for genomic data',
        epilog='Use "aindex <command> --help" for command-specific help'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Add subcommands
    subparsers.add_parser('compute-aindex', help='Compute aindex for k-mer analysis')
    subparsers.add_parser('compute-index', help='Compute index from input data') 
    subparsers.add_parser('compute-reads', help='Process reads')
    subparsers.add_parser('count', help='Count k-mers using fast built-in counter')
    subparsers.add_parser('build-hash', help='Build perfect hash for k-mers')
    subparsers.add_parser('generate', help='Generate all possible k-mers')
    subparsers.add_parser('reads-to-fasta', help='Convert reads to FASTA format')
    subparsers.add_parser('version', help='Show version information')
    subparsers.add_parser('info', help='Show system and installation information')
    
    # Parse main args
    if len(sys.argv) == 1:
        parser.print_help()
        return 1
    
    # Handle special case where user wants help for a subcommand
    if len(sys.argv) >= 3 and sys.argv[2] in ['-h', '--help']:
        # Pass help to subcommand
        args, remaining = parser.parse_known_args(sys.argv[1:2])  # Only parse the command
        remaining = sys.argv[2:]  # Include the --help
    else:
        args, remaining = parser.parse_known_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # Route to appropriate command handler
    command_map = {
        'compute-aindex': cmd_compute_aindex,
        'compute-index': cmd_compute_index,
        'compute-reads': cmd_compute_reads,
        'count': cmd_count_kmers,
        'build-hash': cmd_build_hash,
        'generate': cmd_generate_kmers,
        'reads-to-fasta': cmd_reads_to_fasta,
        'version': cmd_version,
        'info': cmd_info,
    }
    
    if args.command in command_map:
        return command_map[args.command](remaining)
    else:
        print(f"Unknown command: {args.command}")
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())
