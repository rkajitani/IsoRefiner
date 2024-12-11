#!/usr/bin/python

import logging
import os
import re
import sys
import polars as pl
from isorefiner.common import func_with_log, run_command, filter_bam

logger = logging.getLogger(__name__)

        
def main(args):
    try:
        # Set variables
        raw_input_gtf = os.path.abspath(args.input_gtf)
        raw_genome_file = os.path.abspath(args.genome)
        raw_reads_files = [os.path.abspath(_) for _ in args.reads]
        out_gtf = os.path.abspath(args.out_gtf)
        n_thread = args.threads
        max_indel = args.max_indel
        max_clip = args.max_clip
        min_idt = args.min_idt
        min_cov = args.min_cov * 100
        min_mean_depth = args.min_mean_depth

        # Create and move to working directory
        work_dir = args.work_dir
        os.makedirs(work_dir, exist_ok=True)
        os.chdir(work_dir)

        # Set logger
        logging.basicConfig(
            filename="log.txt",
            level=logging.INFO,
            format="%(asctime)s - PID=%(process)d - %(levelname)s - %(message)s"
        )

        # Set input files and create symbolic links
        input_gtf = "input.gtf"
        if os.path.lexists(input_gtf):
            input_gtf = f"raw_{os.getpid()}.gtf"

        genome_file = "genome.fasta"
        if os.path.lexists(genome_file):
            genome_file = f"genome_{os.getpid()}.fasta"
        os.symlink(raw_genome_file, genome_file)

        reads_files = list()
        for i, raw_reads_file in enumerate(raw_reads_files, start=1):
            reads_ext = os.path.basename(raw_reads_file).split(".", 1)[1]
            reads_file = f"reads_{i}.{reads_ext}"
            if os.path.lexists(reads_file):
                reads_file = f"reads_{i}_{os.getpid()}.{reads_ext}"
            reads_files.append(reads_file)
            os.symlink(raw_reads_file, reads_file)

        # Main process
        logger.info(f"Starting isorefiner filter")
        fix_gtf_strand(raw_input_gtf, input_gtf)
        run_command(f"gffread -w asm.fa -g {genome_file} {input_gtf}")
        run_command(f"minimap2 -ax map-ont --secondary=no -t {n_thread} asm.fa {' '.join(reads_files)} | samtools view -b -F 2308 -", stdout="raw.bam")
        run_command(f"samtools sort -@ {n_thread} -m 2G raw.bam", stdout="sorted.bam", stderr="samtools_sort.stderr")
        filter_bam("sorted.bam", "filt.bam", max_indel, max_clip, min_idt)
        run_command(f"samtools index filt.bam")
        run_command(f"samtools coverage filt.bam", stdout="cov.tsv")
        filter_coverage_table("cov.tsv", "cov_depth_filt.list", min_cov, min_mean_depth)
        gtf_extract_transcript(input_gtf, out_gtf, "cov_depth_filt.list")
        logger.info(f"Finished isorefiner filter")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)


@func_with_log
def fix_gtf_strand(in_gtf, out_gtf):
    fix_flag = False
    with open(in_gtf) as fin:
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if f[6] != "+" and f[6] != "-":
                fix_flag = True
                break

    if fix_flag:
        with open(in_gtf) as fin, open(out_gtf, "w") as fout:
            for ln in fin:
                if len(ln) == 0 or ln[0] == "#":
                    continue
                f = ln.rstrip("\n").split("\t")
                if f[6] != "+" and f[6] != "-":
                    f[6] = "+"
                print("\t".join(f), file=fout)
    else:
        os.symlink(in_gtf, out_gtf)


@func_with_log
def filter_coverage_table(in_file, out_file, min_cov, min_mean_depth):
    df = (pl.read_csv(in_file, separator="\t")
        .select("#rname", "coverage", "meandepth")
        .filter((pl.col("coverage") >= min_cov) & (pl.col("meandepth") >= min_mean_depth))
        .select("#rname")
    )
    df.write_csv(out_file, include_header=False)


@func_with_log
def gtf_extract_transcript(in_gtf, out_gtf, target_list_file):
    with open(target_list_file) as fin:
        transcript_set = set(fin.read().splitlines())

    gene_set = set()
    gene_re = re.compile(r'gene_id "([^"]*)"')
    transcript_re = re.compile(r'transcript_id "([^"]*)"')
    with open(in_gtf) as fin:
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if f[2] != "transcript":
                continue
            m = gene_re.search(f[8])
            if not m:
                continue
            gene_id = m.group(1)
            m = transcript_re.search(f[8])
            if not m:
                continue
            transcript_id = m.group(1)
            if transcript_id in transcript_set:
                gene_set.add(gene_id)

    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                print(ln, end="", file=fout)
                continue
            f = ln.rstrip("\n").split("\t")
            if f[2] == "gene":
                m = gene_re.search(f[8])
                if not m:
                    continue
                gene_id = m.group(1)
                if gene_id in gene_set:
                    print(ln, end="", file=fout)
            else:
                m = transcript_re.search(f[8])
                if not m:
                    continue
                transcript_id = m.group(1)
                if transcript_id in transcript_set:
                    print(ln, end="", file=fout)
