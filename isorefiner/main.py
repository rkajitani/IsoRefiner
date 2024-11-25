#!/usr/bin/env python

import argparse
from isorefiner import filter, refine 


def main():
    parser = argparse.ArgumentParser(
        description="IsoRefiner, a tool to filter, merge, and refine transcript isoform structures.",
        formatter_class=argparse.MetavarTypeHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Filter subcommand
    parser_filter = subparsers.add_parser('filter', help='Filter transcript isoforms.')
    parser_filter.add_argument("input.gtf", type=str, help="Input transcript isoform structures (GTF)")
    parser_filter.add_argument("reads.fastq/a(.gz)", type=str, help="Reads (FASTA or FASTQ, gzip allowed)")
    parser_filter.add_argument("genome.fasta", type=str, help="Reference genome (FASTA)")
    parser_filter.add_argument('--max_indel', type=int, default=20, help='Max indel for read mapping')
    parser_filter.add_argument('--max_clip', type=int, default=200, help='Max clip (unaligned) length for read mapping')
    parser_filter.add_argument('--min_idt', type=float, default=0.90, help='Min identity for read mapping [0-1]')
    parser_filter.add_argument('--min_cov', type=float, default=0.95, help='Min coverage for filtering [0-1]')
    parser_filter.add_argument('--min_mean_depth', type=float, default=1.0, help='Min mean coverage depth for filtering')
    parser_filter.set_defaults(func=filter.main)

    # Refine subcomand
    parser_refine = subparsers.add_parser('refine', help='Merge and refine transcript isoforms.')
    parser_refine.set_defaults(func=refine.main)

    args = parser.parse_args()
    args.func(args)
