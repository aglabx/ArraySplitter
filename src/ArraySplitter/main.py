#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ArraySplitter - Main entry point with subcommands.
"""

import argparse
import sys


def run_it():
    """Main entry point with subcommands."""
    
    parser = argparse.ArgumentParser(
        description="ArraySplitter - De novo decomposition and analysis of satellite DNA arrays",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Split arrays into monomers (default action)
  arraysplitter split -i arrays.fa -o output_prefix
  
  # Split with predefined cuts
  arraysplitter split -i arrays.fa -o output_prefix -c ATG,CGCG
  
  # Classify arrays into families (uses .lengths file from split step)
  arraysplitter classify -i output_prefix.lengths -o classification
  
  # Rotate monomers
  arraysplitter rotate -i decomposed.fa -o rotated.fa
  
  # Extract unique monomers
  arraysplitter extract -i decomposed.fa -o monomers

  # Refine decomposition (split merged, join fragments)
  arraysplitter refine -i decomposed.fa -o refined.fa

  # Refine MSA by removing sequences with rare indels
  arraysplitter msa_refine -i aligned.fa -o cleaned

For help on specific commands:
  arraysplitter <command> -h
"""
    )
    
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands'
    )
    
    # Split command (decomposition)
    parser_split = subparsers.add_parser(
        'split',
        help='Split arrays into monomers (decomposition)'
    )
    parser_split.add_argument("-i", "--input", help="Input file", required=True)
    parser_split.add_argument("-o", "--output", help="Output prefix", required=True)
    parser_split.add_argument(
        "--format",
        help="Input format: fasta, trf, satellome [fasta]",
        default="fasta",
        choices=["fasta", "trf", "satellome"]
    )
    parser_split.add_argument(
        "-c", "--cuts",
        help="Comma-separated list of predefined cut sequences (e.g., ATG,ATGATG)",
        default=None
    )
    parser_split.add_argument(
        "-d", "--depth",
        help="Depth for hint discovery (default: 100)",
        type=int,
        default=100
    )
    parser_split.add_argument(
        "-t", "--threads",
        help="Number of threads (currently not used)",
        type=int,
        default=4
    )
    parser_split.add_argument(
        "-v", "--verbose",
        help="Verbose output",
        action="store_true"
    )
    
    # Classify command
    parser_classify = subparsers.add_parser(
        'classify',
        help='Classify arrays into families based on decomposition patterns'
    )
    parser_classify.add_argument("-i", "--input", help="Input .lengths file from decomposition", required=True)
    parser_classify.add_argument("-o", "--output", help="Output prefix", required=True)
    parser_classify.add_argument(
        "-s", "--similarity",
        help="Similarity threshold for clustering (0-1, default: 0.8)",
        type=float,
        default=0.8
    )
    parser_classify.add_argument(
        "-v", "--verbose",
        help="Verbose output",
        action="store_true"
    )
    
    # Rotate command
    parser_rotate = subparsers.add_parser(
        'rotate',
        help='Rotate monomers to start from the same position'
    )
    parser_rotate.add_argument("-i", "--fasta", help="Input FASTA file", required=True)
    parser_rotate.add_argument("-o", "--output", help="Output file", required=True)
    parser_rotate.add_argument(
        "-s", "--start",
        help="Starting kmer for rotation [None]",
        default=None
    )
    
    # Extract command
    parser_extract = subparsers.add_parser(
        'extract',
        help='Extract and count unique monomers'
    )
    parser_extract.add_argument("-i", "--fasta", help="Input FASTA file", required=True)
    parser_extract.add_argument("-o", "--output", help="Output prefix", required=True)
    parser_extract.add_argument(
        "-m", "--min",
        help="Minimum monomer frequency [1]",
        type=int,
        default=1
    )

    # Refine command
    parser_refine = subparsers.add_parser(
        'refine',
        help='Refine decomposition: split merged monomers, join fragments'
    )
    parser_refine.add_argument("-i", "--input", help="Input .monomers.tsv file", required=True)
    parser_refine.add_argument("-o", "--output", help="Output .monomers.tsv file (refined)", required=True)
    parser_refine.add_argument(
        "-l", "--length",
        help="Expected monomer length (auto-detected if not specified)",
        type=int,
        default=None
    )
    parser_refine.add_argument(
        "-m", "--mismatches",
        help="Maximum mismatches for cut sequence matching (default: 3)",
        type=int,
        default=3
    )
    parser_refine.add_argument(
        "-v", "--verbose",
        help="Verbose output",
        action="store_true"
    )

    # MSA Refine command
    parser_msa_refine = subparsers.add_parser(
        'msa_refine',
        help='Remove sequences with rare indels from MSA'
    )
    parser_msa_refine.add_argument("-i", "--input", help="Input aligned FASTA file", required=True)
    parser_msa_refine.add_argument("-o", "--output", help="Output prefix", required=True)
    parser_msa_refine.add_argument(
        "-f", "--min-freq",
        help="Minimum indel frequency threshold (0-1, default: 0.05 = 5%%)",
        type=float,
        default=0.05
    )
    parser_msa_refine.add_argument(
        "-v", "--verbose",
        help="Verbose output",
        action="store_true"
    )

    # Parse arguments
    args = parser.parse_args()
    
    # Show help if no command specified
    if not args.command:
        parser.print_help()
        sys.exit(0)
    
    # Import and run the appropriate module
    if args.command == 'split':
        from .decompose import main as decompose_main
        import os
        
        if not os.path.isfile(args.input):
            print(f"Error: Input file {args.input} not found")
            sys.exit(1)
        
        # Parse cuts if provided
        predefined_cuts = args.cuts.split(',') if args.cuts else None
        
        decompose_main(
            args.input,
            args.output,
            args.format,
            args.threads,
            predefined_cuts,
            args.depth,
            args.verbose
        )
    
    elif args.command == 'classify':
        from .classify import classify_arrays
        import os
        
        if not os.path.isfile(args.input):
            print(f"Error: Input file {args.input} not found")
            sys.exit(1)
        
        if not args.input.endswith('.lengths'):
            print(f"Warning: Input file should be a .lengths file from ArraySplitter decomposition")
        
        classify_arrays(
            args.input,
            args.output,
            similarity_threshold=args.similarity,
            verbose=args.verbose
        )
    
    elif args.command == 'rotate':
        from .rotate import main as rotate_main
        
        rotate_main(args.fasta, args.output, starting_kmer=args.start)
    
    elif args.command == 'extract':
        from .extract import main as extract_main

        extract_main(args.fasta, args.output, min_tf=args.min)

    elif args.command == 'refine':
        from .refine import main as refine_main
        import os

        if not os.path.isfile(args.input):
            print(f"Error: Input file {args.input} not found")
            sys.exit(1)

        refine_main(
            input_tsv=args.input,
            output_tsv=args.output,
            expected_length=args.length,
            max_mismatches=args.mismatches,
            verbose=args.verbose
        )

    elif args.command == 'msa_refine':
        from .msa_refine import refine_msa
        import os

        if not os.path.isfile(args.input):
            print(f"Error: Input file {args.input} not found")
            sys.exit(1)

        output_clean = f"{args.output}.clean.fasta"
        output_removed = f"{args.output}.removed.fasta"
        output_stats = f"{args.output}.indel_stats.tsv"

        stats = refine_msa(
            input_fasta=args.input,
            output_clean=output_clean,
            output_removed=output_removed,
            output_stats=output_stats,
            min_freq=args.min_freq,
            verbose=args.verbose
        )

        print(f"\nRefinement complete:")
        print(f"  Total: {stats['total']} -> Clean: {stats['clean']}, Removed: {stats['removed']}")
        print(f"  Indels: {stats['total_indels']} total, {stats['rare_indels']} rare (< {stats['min_freq']*100:.1f}%)")

    else:
        print(f"Error: Unknown command '{args.command}'")
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    run_it()