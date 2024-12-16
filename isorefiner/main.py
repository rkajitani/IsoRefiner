#!/usr/bin/env python

import argparse
from isorefiner import filter, refine, trim, map, run_stringtie


def main():
    parser = argparse.ArgumentParser(description="IsoRefiner, a tool to filter, merge, and refine transcript isoform structures.")
    subparsers = parser.add_subparsers(dest='command', required=True)

    # Trim subcomand
    parser_trim = subparsers.add_parser("trim", help="Trim nanopore reads using Porechop_ABI.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_trim.add_argument("-r", "--reads", type=str, required=True, nargs="*", help="Reads (FASTQ or FASTA, gzip allowed, mandatory)")
    parser_trim.add_argument("-o", "--out_prefix", type=str, default="isorefiner_trimmed", help="Prefix of final output files (extentions are those of input files)")
    parser_trim.add_argument("-d", "--work_dir", type=str, default="isorefiner_trim_work", help="Working directory containing intermediate and log files")
    parser_trim.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    parser_trim.add_argument("-p", "--tool_option", type=str, default="", help="Option for Porechomp_ABI (quoted string)")
    parser_trim.set_defaults(func=trim.main)

    # Map subcomand
    parser_map = subparsers.add_parser("map", help="Map reads to the reference genome using Minimap2, and sort BAM files.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_map.add_argument("-r", "--reads", type=str, required=True, nargs="*", help="Reads (FASTQ or FASTA, gzip allowed, mandatory)")
    parser_map.add_argument("-g", "--genome", type=str, required=True, help="Reference genome (FASTA, mandatory)")
    parser_map.add_argument("-o", "--out_prefix", type=str, default="isorefiner_mapped", help="Prefix of output BAM files")
    parser_map.add_argument("-d", "--work_dir", type=str, default="isorefiner_map_work", help="Working directory containing intermediate and log files")
    parser_map.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    parser_map.add_argument("-m", "--mm2_option", type=str, default="-x splice -ub -k14 --secondary=no", help="Option for minimap2 (quoted string)")
    parser_map.add_argument("-s", "--sort_option", type=str, default="-m 2G", help="Option for samtools sort (quoted string)")
    parser_map.set_defaults(func=map.main)

    # StringTie subcomand
    parser_run_stringtie = subparsers.add_parser("run_stringtie", help="Run StringTie (read mapping-based tool).", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_run_stringtie.add_argument("-b", "--bam", type=str, required=True, nargs="*", help="Mapped reads files (BAM, mandatory)")
    parser_run_stringtie.add_argument("-g", "--genome", type=str, required=True, help="Reference genome (FASTA, mandatory)")
    parser_run_stringtie.add_argument("-a", "--ref_gtf", type=str, required=True, help="Reference genome annotation (GTF, mandatory)")
    parser_run_stringtie.add_argument("-o", "--out_gtf", type=str, default="isorefiner_stringtie.gtf", help="Final output file name (GTF)")
    parser_run_stringtie.add_argument("-d", "--work_dir", type=str, default="isorefiner_stringtie_work", help="Working directory containing intermediate and log files")
    parser_run_stringtie.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    parser_run_stringtie.add_argument("-p", "--tool_option", type=str, default="", help="Option for StringTie (quoted string)")
    parser_run_stringtie.set_defaults(func=run_stringtie.main)
    
    # Filter subcommand
    parser_filter = subparsers.add_parser("filter", help="Filter transcript isoforms.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_filter.add_argument("-i", "--input_gtf", type=str, required=True, help="Input transcript isoform structures (GTF, mandatory)")
    parser_filter.add_argument("-r", "--reads", type=str, required=True, nargs="*", help="Reads (FASTQ or FASTA, gzip allowed, mandatory)")
    parser_filter.add_argument("-g", "--genome", type=str, required=True, help="Reference genome (FASTA, mandatory)")
    parser_filter.add_argument("-o", "--out_gtf", type=str, default="isorefiner_filtered.gtf", help="Final output file name (GTF)")
    parser_filter.add_argument("-d", "--work_dir", type=str, default="isorefiner_filter_work", help="Working directory containing intermediate and log files")
    parser_filter.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    parser_filter.add_argument("--max_indel", type=int, default=20, help="Max indel for read mapping")
    parser_filter.add_argument("--max_clip", type=int, default=200, help="Max clip (unaligned) length for read mapping")
    parser_filter.add_argument("--min_idt", type=float, default=0.90, help="Min identity for read mapping [0-1]")
    parser_filter.add_argument("--min_cov", type=float, default=0.95, help="Min coverage for filtering [0-1]")
    parser_filter.add_argument("--min_mean_depth", type=float, default=1.0, help="Min mean coverage depth for filtering")
    parser_filter.set_defaults(func=filter.main)

    # Refine subcomand
    parser_refine = subparsers.add_parser("refine", help="Merge and refine transcript isoforms.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_refine.add_argument("-i", "--input_gtf", type=str, required=True, nargs="*", help="Input transcript isoform structures (GTF, mandatory)")
    parser_refine.add_argument("-r", "--reads", type=str, required=True, nargs="*", help="Reads (FASTQ or FASTA, gzip allowed, mandatory)")
    parser_refine.add_argument("-g", "--genome", type=str, required=True, help="Reference genome (FASTA, mandatory)")
    parser_refine.add_argument("-a", "--ref_gtf", type=str, required=True, help="Reference genome annotation (GTF, mandatory)")
    parser_refine.add_argument("-o", "--out_gtf", type=str, default="isorefiner_refined.gtf", help="Final output file name (GTF)")
    parser_refine.add_argument("-d", "--work_dir", type=str, default="isorefiner_refine_work", help="Working directory containing intermediate and log files")
    parser_refine.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    parser_refine.add_argument("--max_indel", type=int, default=20, help="Max indel for read mapping")
    parser_refine.add_argument("--max_clip", type=int, default=200, help="Max clip (unaligned) length for read mapping")
    parser_refine.add_argument("--min_idt", type=float, default=0.90, help="Min identity for read mapping [0-1]")
    parser_refine.add_argument("--intron_dist_th", type=int, default=20, help="Intron distance threshold to exclude erroneous isoforms")
    parser_refine.set_defaults(func=refine.main)

    args = parser.parse_args()
    args.func(args)
