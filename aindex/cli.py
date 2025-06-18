#!/usr/bin/env python3
"""
aindex command-line interface
Unified CLI for all aindex tools and utilities
"""

import argparse
import sys
import os
import platform
import subprocess
import shutil
from pathlib import Path
import importlib.util

# Use importlib instead of deprecated pkg_resources
try:
    import importlib.metadata as importlib_metadata
except ImportError:
    import importlib_metadata


def validate_input_output_files(input_file, output_file, command_name):
    """Validate that input and output files are different to prevent data loss"""
    if not input_file or not output_file:
        return True  # Skip validation if either is empty
    
    try:
        # Convert to Path objects for robust comparison
        input_path = Path(input_file).resolve()
        output_path = Path(output_file).resolve()
        
        # Check if they are the same file
        if input_path == output_path:
            print(f"Error in {command_name}: Input and output files cannot be the same!")
            print(f"Input file:  {input_file}")
            print(f"Output file: {output_file}")
            print("This would overwrite your input data. Please specify a different output file.")
            return False
            
        # Check if output file is a prefix of input file (e.g., input.fastq vs input.fastq.counts)
        # This is usually okay, but warn the user
        if str(input_path).startswith(str(output_path)) or str(output_path).startswith(str(input_path)):
            print(f"Warning in {command_name}: Input and output files have similar names:")
            print(f"Input file:  {input_file}")
            print(f"Output file: {output_file}")
            print("Please verify this is intentional to avoid confusion.")
            
    except Exception as e:
        # If we can't resolve paths, just warn and continue
        print(f"Warning: Could not validate input/output paths: {e}")
    
    return True


def validate_output_file_overwrite(output_file, command_name, force_overwrite=False):
    """Check if output file exists and warn about overwriting"""
    if not output_file:
        return True
    
    try:
        output_path = Path(output_file)
        if output_path.exists():
            if force_overwrite:
                print(f"Warning in {command_name}: Overwriting existing file {output_file}")
                return True
            else:
                print(f"Warning in {command_name}: Output file already exists: {output_file}")
                print("This will overwrite the existing file. Use --force to suppress this warning.")
                # For now, we'll continue but warn the user
                # In the future, we could add --force flag to each command
                return True
    except Exception as e:
        print(f"Warning: Could not check output file: {e}")
    
    return True


def detect_platform():
    """Detect current platform and return optimization info"""
    system = platform.system()
    machine = platform.machine()
    
    platform_info = {
        'system': system,
        'machine': machine,
        'is_apple_silicon': system == 'Darwin' and machine == 'arm64',
        'is_macos': system == 'Darwin',
        'is_linux': system == 'Linux',
        'is_windows': system == 'Windows',
        'cpu_count': os.cpu_count() or 4
    }
    
    return platform_info


def get_optimal_executable(base_name, platform_info=None):
    """Get the optimal executable name based on platform"""
    if platform_info is None:
        platform_info = detect_platform()
    
    # Define platform-specific executable mappings
    # Note: On ARM64, the main binaries are already ARM64-optimized
    executable_variants = {
        'kmer_counter': {
            'default': 'kmer_counter'  # Now always optimal for current platform
        },
        'count_kmers': {
            'default': 'kmer_counter'  # Use the same binary
        },
        # Add more mappings as needed
        'compute_index': {
            'default': 'compute_index'
        },
        'compute_aindex': {
            'default': 'compute_aindex'
        },
        'compute_reads': {
            'default': 'compute_reads'
        },
        'generate_all_13mers': {
            'default': 'generate_all_13mers'
        },
        'build_13mer_hash': {
            'default': 'build_13mer_hash'
        },
        'count_kmers13': {
            'default': 'count_kmers13'
        },
        'compute_aindex13': {
            'default': 'compute_aindex13'
        }
    }
    
    if base_name not in executable_variants:
        return base_name
    
    variants = executable_variants[base_name]
    
    # Always use default since binaries are now platform-optimized at build time
    return variants['default']


def print_platform_info():
    """Print current platform information"""
    platform_info = detect_platform()
    print(f"Platform: {platform_info['system']} {platform_info['machine']}")
    print(f"CPU cores: {platform_info['cpu_count']}")
    
    if platform_info['is_apple_silicon']:
        print("✓ Apple Silicon (ARM64) optimizations available")
    elif platform_info['is_macos']:
        print("• macOS (Intel) - standard optimizations")
    elif platform_info['is_linux']:
        print("• Linux - standard optimizations")
    elif platform_info['is_windows']:
        print("• Windows - standard optimizations")
    
    print()


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
    
    # Method 2: Try via importlib (installed package)
    try:
        try:
            # Python 3.8+
            from importlib.resources import files
            package_path = files('aindex') / 'bin'
        except ImportError:
            # Fallback for older Python versions
            import importlib_resources
            package_path = importlib_resources.files('aindex') / 'bin'
        
        if package_path.exists() and any(package_path.iterdir()):
            return package_path
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


def run_executable(exe_name, args, verbose=False):
    """Run a binary executable from bin directory with platform optimization"""
    platform_info = detect_platform()
    bin_dir = get_bin_path()
    
    # Get optimal executable for current platform
    optimal_exe = get_optimal_executable(exe_name, platform_info)
    
    if verbose:
        print(f"Platform: {platform_info['system']} {platform_info['machine']}")
        print(f"Selected executable: {optimal_exe}")
    
    # Try different executable extensions based on platform
    candidates = [optimal_exe]
    if optimal_exe.endswith('.exe'):
        # If requested with .exe, also try without extension
        candidates.append(optimal_exe[:-4])
    else:
        # If requested without extension, also try with .exe on Windows
        if platform_info['is_windows']:
            candidates.append(optimal_exe + '.exe')
    
    # Also try the original exe_name as fallback
    if exe_name != optimal_exe:
        candidates.append(exe_name)
        if platform_info['is_windows'] and not exe_name.endswith('.exe'):
            candidates.append(exe_name + '.exe')
    
    exe_path = None
    for candidate in candidates:
        candidate_path = bin_dir / candidate
        if candidate_path.exists():
            exe_path = candidate_path
            if verbose:
                print(f"Found executable: {exe_path}")
            break
    
    if exe_path is None:
        print(f"Error: Executable {exe_name} not found in {bin_dir}")
        print(f"Platform: {platform_info['system']} {platform_info['machine']}")
        if platform_info['is_apple_silicon']:
            print("Note: ARM64-optimized version expected but not found")
        print(f"Tried: {[str(bin_dir / c) for c in candidates]}")
        return 1
    
    try:
        # Add platform-specific optimizations to args if needed
        optimized_args = args.copy()
        
        # For ARM64 version, ensure we use all available cores by default
        if 'arm64' in exe_path.name and platform_info['is_apple_silicon']:
            # Check if threads argument is missing and executable supports it
            if exe_path.name in ['kmer_counter_arm64'] and '-t' not in args:
                # Add optimal thread count for M1/M2
                optimized_args.extend(['-t', str(min(platform_info['cpu_count'], 8))])
        
        if verbose:
            print(f"Running: {exe_path.name} {' '.join(optimized_args)}")
        
        # Run the executable with the provided arguments
        result = subprocess.run([str(exe_path)] + optimized_args, check=False)
        return result.returncode
    except Exception as e:
        print(f"Error running {exe_path.name}: {e}")
        return 1


def run_python_script(script_name, args):
    """Run a Python script from scripts or bin directory"""
    # Try to find the script in multiple locations
    current_dir = Path(__file__).parent.parent
    search_dirs = [
        current_dir / "scripts",  # Development environment
        get_bin_path(),           # Installed package bin
        current_dir / "bin",      # Local bin directory
    ]
    
    script_path = None
    for search_dir in search_dirs:
        if search_dir.exists():
            candidate_path = search_dir / script_name
            if candidate_path.exists():
                script_path = candidate_path
                break
    
    if script_path is None:
        print(f"Error: Script {script_name} not found")
        print(f"Searched in: {[str(d) for d in search_dirs if d.exists()]}")
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
    
    # Validate input and output files are different
    if not validate_input_output_files(parsed_args.input, parsed_args.output, 'compute-aindex'):
        return 1
    
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
        description='Compute LU index for reads with perfect hash'
    )
    parser.add_argument('dat_file', help='Data file (or use "dummy" with --mock)')
    parser.add_argument('hash_file', help='Perfect hash file (.hash)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--mock', action='store_true', help='Mock data file flag (use 1)')
    
    parsed_args = parser.parse_args(args)
    
    # compute_index expects: dat_file pf_file output_prefix nthreads mock_flag
    exe_args = [
        parsed_args.dat_file,
        parsed_args.hash_file,
        parsed_args.output,
        str(parsed_args.threads),
        "1" if parsed_args.mock else "0"
    ]
    return run_executable('compute_index', exe_args)


def detect_file_format(filepath):
    """Auto-detect file format by reading the first line"""
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
            
        if not first_line:
            return 'unknown', f"Empty file: {filepath}"
            
        if first_line.startswith('>'):
            return 'fasta', first_line
        elif first_line.startswith('@'):
            return 'fastq', first_line  
        elif all(c in 'ATCGN' for c in first_line.upper()):
            return 'reads', first_line
        else:
            return 'unknown', first_line
            
    except Exception as e:
        return 'error', str(e)


def cmd_compute_reads(args):
    """Process reads using compute_reads"""
    parser = argparse.ArgumentParser(
        prog='aindex compute-reads',
        description='Convert FASTA or FASTQ reads to simple reads format'
    )
    
    # Input options - mutually exclusive groups
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', help='Single input file (FASTA/FASTQ/reads/single-end FASTQ)')
    input_group.add_argument('-1', '--input1', help='First file for paired-end FASTQ')
    
    parser.add_argument('-2', '--input2', help='Second file for paired-end FASTQ (required with -1)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('--format', choices=['fastq', 'fasta', 'se', 'reads'], 
                       help='Force file format (auto-detected if not specified)')
    
    parsed_args = parser.parse_args(args)
    
    # Validate arguments
    if parsed_args.input1 and not parsed_args.input2:
        print("Error: -2/--input2 is required when using -1/--input1")
        return 1
    
    if parsed_args.input2 and not parsed_args.input1:
        print("Error: -1/--input1 is required when using -2/--input2")
        return 1
    
    # Validate input and output files are different
    input_file = parsed_args.input or parsed_args.input1
    if input_file and not validate_input_output_files(input_file, parsed_args.output, 'compute-reads'):
        return 1
    
    # Also check second input file if paired-end
    if parsed_args.input2 and not validate_input_output_files(parsed_args.input2, parsed_args.output, 'compute-reads'):
        return 1
    
    # Determine mode and files
    if parsed_args.input:
        # Single file mode
        file1 = parsed_args.input
        file2 = "-"
        is_paired = False
    else:
        # Paired-end mode
        file1 = parsed_args.input1
        file2 = parsed_args.input2
        is_paired = True
    
    # Auto-detect format if not specified
    if parsed_args.format:
        file_format = parsed_args.format
        print(f"Using specified format: {file_format}")
    else:
        # Auto-detect format from first file
        detected_format, first_line = detect_file_format(file1)
        
        if detected_format == 'error':
            print(f"Error reading file {file1}: {first_line}")
            return 1
        elif detected_format == 'unknown':
            print(f"Unknown format for file {file1}")
            print(f"First line: {first_line}")
            return 1
        else:
            file_format = detected_format
            print(f"Auto-detected format: {file_format}")
            print(f"First line: {first_line}")
    
    # Determine mode and format for compute_reads binary
    if file_format == 'reads':
        mode = 'raw reads'
        compute_format = 'reads'
    elif file_format == 'fasta':
        mode = 'FASTA'
        compute_format = 'fasta'
    elif file_format == 'fastq' or file_format == 'se':
        if is_paired:
            mode = 'PE fastq'
            compute_format = 'fastq'
        else:
            mode = 'SE fastq'
            compute_format = 'se'
    else:
        print(f"Unsupported format: {file_format}")
        return 1
    
    print(f"Processing reads:")
    print(f"  Mode: {mode}")
    print(f"  File 1: {file1}")
    print(f"  File 2: {file2}")
    print(f"  Format: {compute_format}")
    print(f"  Output: {parsed_args.output}")
    
    # compute_reads expects: fastq_file1|fasta_file1|reads_file fastq_file2|- fastq|fasta|se|reads output_prefix
    exe_args = [file1, file2, compute_format, parsed_args.output]
    return run_executable('compute_reads', exe_args)


def cmd_count_kmers(args):
    """Count k-mers using fast built-in counter with platform optimization"""
    parser = argparse.ArgumentParser(
        prog='aindex count',
        description='Count k-mers using built-in fast counter (auto-optimized for current platform)'
    )
    parser.add_argument('-i', '--input', required=True, help='Input FASTA/FASTQ file')
    parser.add_argument('--hash-file', required=True, help='Precomputed perfect hash file (.hash)')
    parser.add_argument('-o', '--output', required=True, help='Output counts file (.tf.bin)')
    parser.add_argument('-k', '--kmer-size', type=int, choices=[13, 23], default=13, help='K-mer size')
    parser.add_argument('-t', '--threads', type=int, default=None, help='Number of threads (default: auto-optimized)')
    parser.add_argument('--verbose', action='store_true', help='Show platform optimization details')
    
    parsed_args = parser.parse_args(args)
    
    # Validate input and output files are different
    if not validate_input_output_files(parsed_args.input, parsed_args.output, 'count'):
        return 1
    
    platform_info = detect_platform()
    if parsed_args.verbose:
        print_platform_info()
    
    # Auto-optimize thread count if not specified
    if parsed_args.threads is None:
        if platform_info['is_apple_silicon']:
            # M1/M2 optimal: use all P-cores + some E-cores
            parsed_args.threads = min(platform_info['cpu_count'], 8)
        else:
            # Standard optimization
            parsed_args.threads = platform_info['cpu_count']
    
    if parsed_args.kmer_size == 13:
        print(f"Counting 13-mers in {parsed_args.input}")
        if platform_info['is_apple_silicon']:
            print("Using ARM64-optimized k-mer counter")
        
        # count_kmers13 expects: input_file hash_file output_tf_file [num_threads]
        exe_args = [parsed_args.input, parsed_args.hash_file, parsed_args.output, str(parsed_args.threads)]
        return run_executable('count_kmers13', exe_args, parsed_args.verbose)
    else:
        print(f"Counting {parsed_args.kmer_size}-mers in {parsed_args.input}")
        if platform_info['is_apple_silicon']:
            print("Using ARM64-optimized k-mer counter")
        
        # For general k-mer counter, use format: input k output [-t threads] [-m min_count]
        exe_args = [parsed_args.input, str(parsed_args.kmer_size), parsed_args.output, 
                   '-t', str(parsed_args.threads)]
        return run_executable('kmer_counter', exe_args, parsed_args.verbose)


def cmd_count_kmers_direct(args):
    """Count k-mers directly from sequences (no hash required) - ARM64 optimized"""
    parser = argparse.ArgumentParser(
        prog='aindex count-direct',
        description='Count k-mers directly from sequences using optimized counter'
    )
    parser.add_argument('-i', '--input', required=True, help='Input FASTA/FASTQ file')
    parser.add_argument('-k', '--kmer-size', type=int, default=13, help='K-mer size')
    parser.add_argument('-o', '--output', required=True, help='Output k-mer counts file')
    parser.add_argument('-t', '--threads', type=int, default=None, help='Number of threads (default: auto-optimized)')
    parser.add_argument('-m', '--min-count', type=int, default=1, help='Minimum count threshold')
    parser.add_argument('--verbose', action='store_true', help='Show platform optimization details')
    
    parsed_args = parser.parse_args(args)
    
    # Validate input and output files are different
    if not validate_input_output_files(parsed_args.input, parsed_args.output, 'count-direct'):
        return 1
    
    platform_info = detect_platform()
    if parsed_args.verbose:
        print_platform_info()
    
    # Auto-optimize thread count if not specified
    if parsed_args.threads is None:
        if platform_info['is_apple_silicon']:
            # M1/M2 optimal: use all cores efficiently
            parsed_args.threads = platform_info['cpu_count']
        else:
            # Standard optimization
            parsed_args.threads = platform_info['cpu_count']
    
    print(f"Direct k-mer counting: k={parsed_args.kmer_size}, input={parsed_args.input}")
    if platform_info['is_apple_silicon']:
        print("Using ARM64-optimized direct k-mer counter")
    
    # Use the optimized k-mer counter directly
    # Format: input k output [-t threads] [-m min_count]
    exe_args = [parsed_args.input, str(parsed_args.kmer_size), parsed_args.output,
               '-t', str(parsed_args.threads), '-m', str(parsed_args.min_count)]
    
    return run_executable('kmer_counter', exe_args, parsed_args.verbose)


def cmd_build_hash(args):
    """Build perfect hash for k-mers"""
    parser = argparse.ArgumentParser(
        prog='aindex build-hash',
        description='Build perfect hash for k-mers'
    )
    parser.add_argument('-i', '--input', required=True, help='Input k-mers file')
    parser.add_argument('-o', '--output', required=True, help='Output hash file')
    parser.add_argument('-k', '--kmer-size', type=int, choices=[13, 23], default=13, help='K-mer size')
    
    parsed_args = parser.parse_args(args)
    
    # Validate input and output files are different
    if not validate_input_output_files(parsed_args.input, parsed_args.output, 'build-hash'):
        return 1
    
    if parsed_args.kmer_size == 13:
        print(f"Building 13-mer hash for {parsed_args.input}")
        # build_13mer_hash expects: kmers_file output_hash_file
        exe_args = [parsed_args.input, parsed_args.output]
        return run_executable('build_13mer_hash', exe_args)
    else:
        print(f"Building hash for {parsed_args.kmer_size}-mers")
        # Use general purpose hash builder
        exe_args = [parsed_args.input, parsed_args.output]
        return run_executable('compute_mphf_seq', exe_args)


def cmd_generate_kmers(args):
    """Generate all possible k-mers"""
    parser = argparse.ArgumentParser(
        prog='aindex generate',
        description='Generate all possible k-mers'
    )
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-i', '--with-indices', action='store_true', help='Include numerical indices in output')
    parser.add_argument('-b', '--binary', action='store_true', help='Generate binary format')
    parser.add_argument('-s', '--stats', action='store_true', help='Show statistics only')
    parser.add_argument('-v', '--validate', action='store_true', help='Run validation test')
    
    parsed_args = parser.parse_args(args)
    
    # Check if output file already exists
    if not validate_output_file_overwrite(parsed_args.output, 'generate'):
        return 1
    
    print(f"Generating all 13-mers to {parsed_args.output}")
    # generate_all_13mers expects: output_file [options]
    exe_args = [parsed_args.output]
    if parsed_args.with_indices:
        exe_args.append('-i')
    if parsed_args.binary:
        exe_args.append('-b')
    if parsed_args.stats:
        exe_args.append('-s')
    if parsed_args.validate:
        exe_args.append('-v')
    return run_executable('generate_all_13mers', exe_args)
   

def cmd_compute_aindex_direct(args):
    """Direct call to compute_aindex binary for expert users"""
    parser = argparse.ArgumentParser(
        prog='aindex compute-aindex-direct',
        description='Direct call to compute_aindex binary (for expert users)'
    )
    parser.add_argument('reads_file', help='Input reads file (one sequence per line)')
    parser.add_argument('hash_file', help='Precomputed perfect hash file (.hash)')
    parser.add_argument('output_prefix', help='Output prefix for generated files')
    parser.add_argument('-t', '--threads', type=int, required=True, help='Number of threads to use')
    parser.add_argument('-k', '--kmer-size', type=int, choices=[13, 23], default=13, help='K-mer size (13 or 23)')
    parser.add_argument('--tf-file', help='TF frequencies file (.tf.bin)')
    parser.add_argument('--kmers-bin', help='Binary k-mers file (.kmers.bin) [for k=23 only]')
    parser.add_argument('--kmers-text', help='Text k-mers file (.kmers) [for k=23 only]')
    
    parsed_args = parser.parse_args(args)
    
    if parsed_args.kmer_size == 13:
        print(f"Computing 13-mer aindex for {parsed_args.reads_file}")
        if not parsed_args.tf_file:
            print("Error: --tf-file is required for 13-mer mode")
            return 1
        # compute_aindex13 expects: reads_file hash_file tf_file output_prefix num_threads
        exe_args = [
            parsed_args.reads_file,
            parsed_args.hash_file, 
            parsed_args.tf_file,
            parsed_args.output_prefix,
            str(parsed_args.threads)
        ]
        return run_executable('compute_aindex13', exe_args)
    else:
        print(f"Computing 23-mer aindex for {parsed_args.reads_file}")
        if not all([parsed_args.tf_file, parsed_args.kmers_bin, parsed_args.kmers_text]):
            print("Error: --tf-file, --kmers-bin, and --kmers-text are required for 23-mer mode")
            return 1
        # compute_aindex expects: reads_file hash_file output_prefix num_threads k tf_file kmers_bin_file kmers_text_file
        exe_args = [
            parsed_args.reads_file,
            parsed_args.hash_file,
            parsed_args.output_prefix,
            str(parsed_args.threads),
            str(parsed_args.kmer_size),
            parsed_args.tf_file,
            parsed_args.kmers_bin,
            parsed_args.kmers_text
        ]
        return run_executable('compute_aindex', exe_args)


def cmd_reads_to_fasta(args):
    """Convert reads to FASTA format"""
    parser = argparse.ArgumentParser(
        prog='aindex reads-to-fasta',
        description='Convert reads to FASTA format'
    )
    parser.add_argument('-i', '--input', required=True, help='Reads input file')
    parser.add_argument('-o', '--output', required=True, help='FASTA output file')
    
    parsed_args = parser.parse_args(args)
    
    # Validate input and output files are different
    if not validate_input_output_files(parsed_args.input, parsed_args.output, 'reads-to-fasta'):
        return 1
    
    print(f"Converting reads from {parsed_args.input} to FASTA format...")
    script_args = ['-i', parsed_args.input, '-o', parsed_args.output]
    return run_python_script('reads_to_fasta.py', script_args)


def cmd_help(args):
    """Show detailed help for all commands"""
    platform_info = detect_platform()
    
    print("=== aindex Command Line Interface ===")
    print()
    print_platform_info()
    print("Available commands:")
    print()
    
    commands = {
        'generate': 'Generate all possible k-mers (13-mers)',
        'build-hash': 'Build perfect hash for k-mers',
        'count': 'Count k-mers using fast built-in counter (with hash)',
        'count-direct': 'Count k-mers directly from sequences (no hash required) ⚡',
        'compute-reads': 'Convert FASTA/FASTQ reads to simple reads format',
        'compute-aindex': 'Compute aindex for k-mer analysis (high-level)',
        'compute-aindex-direct': 'Direct call to compute_aindex binary (expert)',
        'compute-index': 'Compute LU index for reads with perfect hash',
        'reads-to-fasta': 'Convert reads to FASTA format',
        'version': 'Show version information',
        'info': 'Show system and installation information (--skip-cpp-test, --verbose-debug, --file-details available)',
        'platform': 'Show platform optimization information',
        'api-docs': 'Show detailed C++ API documentation'
    }
    
    for cmd, desc in commands.items():
        if cmd == 'count-direct' and platform_info['is_apple_silicon']:
            print(f"  {cmd:<25} {desc} (ARM64-optimized)")
        else:
            print(f"  {cmd:<25} {desc}")
    
    print()
    print("Typical workflows:")
    print()
    print("📈 Fast direct k-mer counting (recommended):")
    print("  aindex count-direct -i input.fastq -k 13 -o kmers.txt --verbose")
    print()
    print("🔬 Traditional workflow with perfect hash (for large datasets):")
    print("  1. aindex generate -o all_13mers.txt")
    print("  2. aindex build-hash -i all_13mers.txt -o 13mer_index.hash")
    print("  3. aindex count -i input.fastq --hash-file 13mer_index.hash -o counts.tf.bin")
    print("  4. aindex compute-reads -i input.fastq -o reads_prefix  # single-end")
    print("     OR: aindex compute-reads -1 R1.fastq -2 R2.fastq -o reads_prefix  # paired-end")
    print("  5. aindex compute-aindex-direct reads_prefix.reads 13mer_index.hash output_prefix -t 4 --tf-file counts.tf.bin")
    print()
    
    if platform_info['is_apple_silicon']:
        print("💡 Apple Silicon optimizations:")
        print("  • ARM64-optimized k-mer counting with NEON instructions")
        print("  • Memory layout optimized for M1/M2 cache hierarchy")
        print("  • Automatic optimal thread count selection")
        print()
    
    print("Use 'aindex <command> --help' for detailed help on specific commands.")
    print("Use 'aindex platform --list-executables' to see all available tools.")
    print("Use 'aindex api-docs' to see detailed C++ API documentation.")
    return 0


def cmd_version(args):
    """Show aindex version"""
    try:
        import aindex
        print(f"aindex version {aindex.__version__}")
        
        # Show available tools
        bin_dir = get_bin_path()
        if bin_dir.exists():
            try:
                executables = []
                for f in bin_dir.iterdir():
                    try:
                        if (f.is_file() and 
                            not f.name.endswith('.py') and 
                            not f.name.startswith('__') and
                            not f.name.startswith('.')):
                            executables.append(f.name)
                    except (OSError, PermissionError):
                        continue  # Skip problematic files silently
                
                if executables:
                    print(f"Available executables in {bin_dir}:")
                    for exe in sorted(executables):
                        print(f"  {exe}")
            except (OSError, PermissionError) as e:
                print(f"Warning: Could not list bin directory: {e}")
        
        return 0
    except ImportError:
        print("aindex not properly installed")
        return 1


def cmd_info(args):
    """Show system and installation information"""
    # Parse args to check for debug options
    parser = argparse.ArgumentParser(
        prog='aindex info',
        description='Show system and installation information'
    )
    parser.add_argument('--skip-cpp-test', action='store_true', 
                       help='Skip C++ API test (use if experiencing crashes)')
    parser.add_argument('--verbose-debug', action='store_true',
                       help='Show very detailed debug information for troubleshooting')
    parser.add_argument('--file-details', action='store_true',
                       help='Show detailed file information (size, permissions, etc.)')
    parser.add_argument('--skip-bin-check', action='store_true',
                       help='Skip bin directory analysis entirely')
    parser.add_argument('--minimal', action='store_true',
                       help='Show only basic information')
    parser.add_argument('--force-gc', action='store_true',
                       help='Force garbage collection to trigger potential issues early')
    parser.add_argument('--test-only', choices=['cpp', 'path', 'bin'], 
                       help='Test only specific component (cpp, path, or bin)')
    
    try:
        parsed_args = parser.parse_args(args)
    except:
        # If parsing fails, proceed with default behavior
        parsed_args = argparse.Namespace(
            skip_cpp_test=False, 
            verbose_debug=False, 
            file_details=False,
            skip_bin_check=False,
            minimal=False,
            force_gc=False,
            test_only=None
        )
    
    try:
        import aindex
        print("=== aindex System Information ===")
        print(f"Version: {aindex.__version__}")
        print(f"Python: {sys.version}")
        print(f"Platform: {sys.platform}")
        
        # Test C++ module safely
        if not parsed_args.test_only or parsed_args.test_only == 'cpp':
            if parsed_args.skip_cpp_test:
                print("C++ API: Test skipped (--skip-cpp-test)")
            else:
                try:
                    print("Testing C++ API...")
                    import aindex.core.aindex_cpp as aindex_cpp
                    wrapper = aindex_cpp.AindexWrapper()
                    methods = [m for m in dir(wrapper) if not m.startswith('_')]
                    print(f"C++ API: Available ({len(methods)} methods)")
                    
                    # Show method descriptions if verbose
                    if parsed_args.verbose_debug:
                        print("\nC++ API Methods:")
                        for method_name in sorted(methods):
                            try:
                                method = getattr(wrapper, method_name)
                                if hasattr(method, '__doc__') and method.__doc__:
                                    doc = method.__doc__.strip()
                                    # Format multi-line docstrings nicely
                                    doc_lines = doc.split('\n')
                                    if len(doc_lines) > 1:
                                        print(f"  {method_name}:")
                                        for line in doc_lines:
                                            if line.strip():
                                                print(f"    {line.strip()}")
                                    else:
                                        print(f"  {method_name}: {doc}")
                                else:
                                    print(f"  {method_name}: (no documentation)")
                            except Exception as e:
                                print(f"  {method_name}: (error getting info: {e})")
                    
                    # Test if we can safely create and destroy multiple instances
                    if parsed_args.verbose_debug:
                        print("\nDEBUG: Testing multiple wrapper instances...")
                        for i in range(3):
                            test_wrapper = aindex_cpp.AindexWrapper()
                            print(f"DEBUG: Created wrapper {i+1}")
                            del test_wrapper
                            print(f"DEBUG: Deleted wrapper {i+1}")
                    
                    # Explicitly delete the wrapper to avoid potential destructor issues
                    print("DEBUG: Deleting main C++ wrapper...")
                    del wrapper
                    if parsed_args.verbose_debug:
                        print("DEBUG: Main wrapper deleted successfully")
                    
                    # Force cleanup of the module reference
                    del aindex_cpp
                    if parsed_args.verbose_debug:
                        print("DEBUG: Module reference deleted")
                        
                    sys.stdout.flush()
                    
                except ImportError as e:
                    print(f"C++ API: Import Error - {e}")
                except Exception as e:
                    print(f"C++ API: Error - {e}")
                    print("Note: Use --skip-cpp-test if this causes crashes")
        
        # Show detailed path information
        if (not parsed_args.minimal and 
            (not parsed_args.test_only or parsed_args.test_only == 'path')):
            print("\n=== Path Information ===")
            
            # Show aindex package location
            try:
                package_dir = Path(aindex.__file__).parent
                print(f"Package location: {package_dir}")
                print(f"Package type: {'Development' if str(package_dir).endswith('workspace/aindex/aindex') else 'Installed'}")
            except Exception as e:
                print(f"Package location: Error - {e}")
        
        # Show bin directory search results with enhanced error handling and detailed debugging
        if (not parsed_args.skip_bin_check and not parsed_args.minimal and 
            (not parsed_args.test_only or parsed_args.test_only == 'bin')):
            try:
                print("DEBUG: Starting bin directory analysis...")
                bin_dir = get_bin_path()
                print(f"Bin directory: {bin_dir}")
                print(f"Bin exists: {bin_dir.exists()}")
                
                if bin_dir.exists():
                    try:
                        print("DEBUG: Listing directory contents...")
                        files = list(bin_dir.iterdir())
                        print(f"Bin files: {len(files)} files")
                        
                        print("DEBUG: Starting file analysis...")
                        # Show executables with safe file checking and detailed debugging
                        executables = []
                        for i, f in enumerate(files):
                            try:
                                if parsed_args.verbose_debug:
                                    print(f"DEBUG: Checking file {i+1}/{len(files)}: {f.name}")
                                
                                # Check if it's a file
                                is_file = f.is_file()
                                if parsed_args.verbose_debug:
                                    print(f"DEBUG:   - is_file(): {is_file}")
                                
                                if parsed_args.file_details:
                                    try:
                                        stat_info = f.stat()
                                        print(f"DEBUG:   - File size: {stat_info.st_size} bytes")
                                        print(f"DEBUG:   - Permissions: {oct(stat_info.st_mode)}")
                                        print(f"DEBUG:   - Is executable: {bool(stat_info.st_mode & 0o111)}")
                                    except Exception as e:
                                        print(f"DEBUG:   - Could not get file stats: {e}")
                                
                                if is_file:
                                    # Check file name conditions
                                    not_py = not f.name.endswith('.py')
                                    not_dunder = not f.name.startswith('__')
                                    not_hidden = not f.name.startswith('.')
                                    
                                    if parsed_args.verbose_debug:
                                        print(f"DEBUG:   - not .py: {not_py}")
                                        print(f"DEBUG:   - not __*: {not_dunder}")
                                        print(f"DEBUG:   - not hidden: {not_hidden}")
                                    
                                    if not_py and not_dunder and not_hidden:
                                        if parsed_args.verbose_debug:
                                            print(f"DEBUG:   - Adding {f.name} to executables")
                                        executables.append(f.name)
                                    else:
                                        if parsed_args.verbose_debug:
                                            print(f"DEBUG:   - Skipping {f.name} (failed conditions)")
                                else:
                                    if parsed_args.verbose_debug:
                                        print(f"DEBUG:   - Skipping {f.name} (not a file)")
                                    
                            except (OSError, PermissionError) as e:
                                print(f"WARNING: Could not check file {f.name}: {e}")
                                continue
                            except Exception as e:
                                print(f"ERROR: Unexpected error checking file {f.name}: {e}")
                                continue
                        
                        if parsed_args.verbose_debug:
                            print("DEBUG: File analysis complete")
                        
                        if executables:
                            if parsed_args.verbose_debug:
                                print("DEBUG: Preparing to display executables list...")
                            sorted_executables = sorted(executables)
                            if parsed_args.verbose_debug:
                                print(f"DEBUG: Sorted {len(sorted_executables)} executables")
                                print("DEBUG: About to join and print executables...")
                                
                            # Try to print executables one by one to isolate the problem
                            if parsed_args.file_details:
                                print("Executables (detailed):")
                                for exe in sorted_executables:
                                    try:
                                        print(f"  - {exe}")
                                        sys.stdout.flush()  # Force flush after each line
                                    except Exception as e:
                                        print(f"ERROR printing executable name {exe}: {e}")
                            else:
                                try:
                                    exe_list = ", ".join(sorted_executables)
                                    if parsed_args.verbose_debug:
                                        print(f"DEBUG: Joined string length: {len(exe_list)}")
                                    print("Executables:", exe_list)
                                    sys.stdout.flush()
                                except Exception as e:
                                    print(f"ERROR creating executables list: {e}")
                                    # Fallback: print one by one
                                    print("Executables (fallback):")
                                    for exe in sorted_executables:
                                        try:
                                            print(f"  - {exe}")
                                            sys.stdout.flush()
                                        except Exception as e:
                                            print(f"ERROR printing {exe}: {e}")
                            
                            if parsed_args.verbose_debug:
                                print("DEBUG: Executables list displayed successfully")
                        else:
                            print("No valid executables found")
                            
                    except (OSError, PermissionError) as e:
                        print(f"Error accessing bin directory: {e}")
                    except Exception as e:
                        print(f"Unexpected error in bin directory processing: {e}")
                else:
                    print("Bin directory does not exist")
                    
                print("DEBUG: Bin directory analysis complete")
            except Exception as e:
                print(f"Error checking bin directory: {e}")
        else:
            if parsed_args.skip_bin_check:
                print("Bin directory check skipped (--skip-bin-check)")
            
        print("DEBUG: About to return from cmd_info function...")
        sys.stdout.flush()  # Force flush output before potential crash
        
        # Explicitly clean up variables that might cause issues on destruction
        try:
            print("DEBUG: Cleaning up variables...")
            if 'bin_dir' in locals():
                del bin_dir
            if 'files' in locals():
                del files
            if 'executables' in locals():
                del executables
            if 'sorted_executables' in locals():
                del sorted_executables
            if 'package_dir' in locals():
                del package_dir
            print("DEBUG: Variables cleaned up")
            sys.stdout.flush()
            
            # Force garbage collection if requested
            if parsed_args.force_gc:
                import gc
                print("DEBUG: Forcing garbage collection...")
                collected = gc.collect()
                print(f"DEBUG: Collected {collected} objects")
                sys.stdout.flush()
                
        except Exception as e:
            print(f"DEBUG: Error during cleanup: {e}")
        
        print("DEBUG: About to return 0...")
        sys.stdout.flush()
        return 0
    except ImportError as e:
        print(f"Error: aindex module not found - {e}")
        return 1
    except Exception as e:
        print(f"Error getting info: {e}")
        return 1


def cmd_api_docs(args):
    """Show detailed C++ API documentation"""
    parser = argparse.ArgumentParser(
        prog='aindex api-docs',
        description='Show detailed C++ API documentation'
    )
    parser.add_argument('--method', type=str, help='Show documentation for specific method')
    parser.add_argument('--category', choices=['loading', 'query', 'utility', 'all'], 
                       default='all', help='Show methods by category')
    parser.add_argument('--examples', action='store_true', help='Show usage examples')
    
    parsed_args = parser.parse_args(args)
    
    try:
        import aindex.core.aindex_cpp as aindex_cpp
        wrapper = aindex_cpp.AindexWrapper()
        
        # Categorize methods
        loading_methods = [
            'load', 'load_hash_file', 'load_reads', 'load_reads_in_memory', 
            'load_aindex', 'load_13mer_index', 'load_13mer_aindex',
            'load_from_prefix_23mer', 'load_aindex_from_prefix_23mer',
            'load_from_prefix_13mer', 'load_aindex_from_prefix_13mer'
        ]
        
        query_methods = [
            'get_tf_value', 'get_tf_values', 'get_hash_value', 'get_hash_values',
            'get_positions', 'get_reads_se_by_kmer', 'get_read_by_rid', 'get_read',
            'get_kid_by_kmer', 'get_kmer_by_kid', 'get_strand', 'get_kmer_info',
            'get_rid', 'get_start', 'get_hash_size', 'get_reads_size',
            'get_tf_value_13mer', 'get_tf_values_13mer', 'get_total_tf_value_13mer',
            'get_tf_both_directions_13mer', 'get_positions_13mer'
        ]
        
        utility_methods = [
            'get_index_info', 'get_13mer_statistics', 'get_23mer_statistics',
            'get_reverse_complement_13mer', 'get_reverse_complement_23mer',
            'debug_kmer_tf_values'
        ]
        
        all_methods = [m for m in dir(wrapper) if not m.startswith('_')]
        
        # Filter methods by category
        if parsed_args.category == 'loading':
            methods_to_show = loading_methods
        elif parsed_args.category == 'query':
            methods_to_show = query_methods
        elif parsed_args.category == 'utility':
            methods_to_show = utility_methods
        else:
            methods_to_show = all_methods
        
        # Filter by specific method if requested
        if parsed_args.method:
            if parsed_args.method in all_methods:
                methods_to_show = [parsed_args.method]
            else:
                print(f"Method '{parsed_args.method}' not found in C++ API")
                print(f"Available methods: {', '.join(sorted(all_methods))}")
                return 1
        
        print("=== aindex C++ API Documentation ===")
        print(f"Total methods: {len(all_methods)}")
        print(f"Showing category: {parsed_args.category}")
        print()
        
        for method_name in sorted(methods_to_show):
            if method_name not in all_methods:
                continue
                
            try:
                method = getattr(wrapper, method_name)
                print(f"📎 {method_name}")
                
                # Show category
                if method_name in loading_methods:
                    print("   Category: Loading")
                elif method_name in query_methods:
                    print("   Category: Query")
                elif method_name in utility_methods:
                    print("   Category: Utility")
                else:
                    print("   Category: Other")
                
                # Show documentation
                if hasattr(method, '__doc__') and method.__doc__:
                    doc = method.__doc__.strip()
                    doc_lines = doc.split('\n')
                    print("   Description:")
                    for line in doc_lines:
                        if line.strip():
                            print(f"     {line.strip()}")
                else:
                    print("   Description: (no documentation available)")
                
                # Show signature if available
                try:
                    import inspect
                    sig = inspect.signature(method)
                    print(f"   Signature: {method_name}{sig}")
                except:
                    pass
                
                # Show usage examples for common methods
                if parsed_args.examples and method_name in ['get_tf_value', 'load_from_prefix_13mer', 'get_positions']:
                    print("   Example:")
                    if method_name == 'get_tf_value':
                        print("     wrapper.get_tf_value('ATCGATCGATCGA')")
                    elif method_name == 'load_from_prefix_13mer':
                        print("     wrapper.load_from_prefix_13mer('data/13mer_index')")
                    elif method_name == 'get_positions':
                        print("     positions = wrapper.get_positions('ATCGATCGATCGA')")
                
                print()
                
            except Exception as e:
                print(f"   Error getting info for {method_name}: {e}")
                print()
        
        # Show category summary
        print("=== Category Summary ===")
        print(f"Loading methods ({len(loading_methods)}): Data loading and initialization")
        print(f"Query methods ({len(query_methods)}): K-mer and read queries") 
        print(f"Utility methods ({len(utility_methods)}): Statistics and utilities")
        print()
        print("Use 'aindex api-docs --category <category>' to filter by category")
        print("Use 'aindex api-docs --method <method_name>' for specific method info")
        print("Use 'aindex api-docs --examples' to see usage examples")
        
        # Clean up
        del wrapper
        del aindex_cpp
        
        return 0
        
    except ImportError as e:
        print(f"Error: Cannot import C++ API - {e}")
        return 1
    except Exception as e:
        print(f"Error getting API documentation: {e}")
        return 1


def cmd_platform_info(args):
    """Show platform information and available optimizations"""
    parser = argparse.ArgumentParser(
        prog='aindex platform',
        description='Show platform information and available optimizations'
    )
    parser.add_argument('--list-executables', action='store_true', 
                       help='List all available executables')
    
    parsed_args = parser.parse_args(args)
    
    platform_info = detect_platform()
    
    print("=== aindex Platform Information ===")
    print_platform_info()
    
    # Show available optimizations
    print("Available optimizations:")
    if platform_info['is_apple_silicon']:
        print("✓ ARM64-optimized binaries built for Apple Silicon")
        print("✓ Apple Silicon memory layout optimizations")
        print("✓ M1/M2 cache-friendly algorithms")
        print("✓ Native ARM64 instruction optimizations")
    else:
        print("• Standard x86_64 optimizations")
        print("• Multi-threading support")
    
    print(f"Recommended thread count: {platform_info['cpu_count']}")
    
    if parsed_args.list_executables:
        print("\n=== Available Executables ===")
        bin_dir = get_bin_path()
        if bin_dir.exists():
            try:
                executables = []
                for file_path in bin_dir.iterdir():
                    try:
                        if (file_path.is_file() and 
                            os.access(file_path, os.R_OK) and  # Check readable instead of executable
                            not file_path.name.endswith('.py') and
                            not file_path.name.startswith('__') and
                            not file_path.name.startswith('.')):
                            executables.append(file_path.name)
                    except (OSError, PermissionError) as e:
                        print(f"Warning: Could not check file {file_path.name}: {e}")
                        continue
                
                if executables:
                    for exe in sorted(executables):
                        if platform_info['is_apple_silicon']:
                            print(f"✓ {exe} (ARM64-optimized)")
                        else:
                            print(f"✓ {exe}")
                else:
                    print("No executables found in bin directory")
            except (OSError, PermissionError) as e:
                print(f"Error accessing bin directory: {e}")
        else:
            print(f"Bin directory not found: {bin_dir}")
    
    return 0


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        prog='aindex',
        description='aindex: perfect hash based index for genomic data',
        epilog='Use "aindex <command> --help" for command-specific help'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Add subcommands
    subparsers.add_parser('help', help='Show detailed help for all commands')
    subparsers.add_parser('generate', help='Generate all possible k-mers (13-mers only)')
    subparsers.add_parser('build-hash', help='Build perfect hash for k-mers')
    subparsers.add_parser('count', help='Count k-mers using fast built-in counter (requires hash file)')
    subparsers.add_parser('count-direct', help='Count k-mers directly from sequences (ARM64-optimized)')
    subparsers.add_parser('compute-reads', help='Convert FASTA/FASTQ reads to simple reads format')
    subparsers.add_parser('compute-aindex', help='Compute aindex for k-mer analysis')
    subparsers.add_parser('compute-aindex-direct', help='Direct call to compute_aindex binary (expert)')
    subparsers.add_parser('compute-index', help='Compute index from input data') 
    subparsers.add_parser('reads-to-fasta', help='Convert reads to FASTA format')
    subparsers.add_parser('version', help='Show version information')
    subparsers.add_parser('info', help='Show system and installation information')
    subparsers.add_parser('platform', help='Show platform information and optimizations')
    subparsers.add_parser('api-docs', help='Show detailed C++ API documentation')
    
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
        'help': cmd_help,
        'generate': cmd_generate_kmers,
        'build-hash': cmd_build_hash,
        'count': cmd_count_kmers,
        'count-direct': cmd_count_kmers_direct,
        'compute-reads': cmd_compute_reads,
        'compute-aindex': cmd_compute_aindex,
        'compute-aindex-direct': cmd_compute_aindex_direct,
        'compute-index': cmd_compute_index,
        'reads-to-fasta': cmd_reads_to_fasta,
        'version': cmd_version,
        'info': cmd_info,
        'platform': cmd_platform_info,
        'api-docs': cmd_api_docs,
    }
    
    if args.command in command_map:
        return command_map[args.command](remaining)
    else:
        print(f"Unknown command: {args.command}")
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())
