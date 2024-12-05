#!/usr/bin/env python

import argparse
from isorefiner import filter, refine 


def main():
    parser = argparse.ArgumentParser(description="IsoRefiner, a tool to filter, merge, and refine transcript isoform structures.")
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Filter subcommand
    parser_filter = subparsers.add_parser("filter", help="Filter transcript isoforms.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_filter.add_argument("input_gtf", type=str, help="Input transcript isoform structures (GTF, mandatory)")
    parser_filter.add_argument("reads_fastq", type=str, help="Reads (FASTQ or FASTA, gzip allowed, mandatory)")
    parser_filter.add_argument("genome_fasta", type=str, help="Reference genome (FASTA, mandatory)")
    parser_filter.add_argument("out_dir", type=str, help="Directory for output and temporary files")
    parser_filter.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    parser_filter.add_argument("--max_indel", type=int, default=20, help="Max indel for read mapping")
    parser_filter.add_argument("--max_clip", type=int, default=200, help="Max clip (unaligned) length for read mapping")
    parser_filter.add_argument("--min_idt", type=float, default=0.90, help="Min identity for read mapping [0-1]")
    parser_filter.add_argument("--min_cov", type=float, default=0.95, help="Min coverage for filtering [0-1]")
    parser_filter.add_argument("--min_mean_depth", type=float, default=1.0, help="Min mean coverage depth for filtering")
    parser_filter.set_defaults(func=filter.main)

    # Refine subcomand
    parser_refine = subparsers.add_parser("refine", help="Merge and refine transcript isoforms.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_refine.set_defaults(func=refine.main)

    args = parser.parse_args()
    args.func(args)
