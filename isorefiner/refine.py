#!/usr/bin/python

import logging
import os
import re
import subprocess
import sys
import polars as pl
from isorefiner.common import func_with_log, run_command, filter_bam

logger = logging.getLogger(__name__)


def main(args):
    try:
        # Set variables
        raw_input_gtfs = [os.path.abspath(_) for _ in args.input_gtf]
        raw_genome_file = os.path.abspath(args.genome)
        raw_ref_gtf = os.path.abspath(args.ref_gtf)
        raw_reads_files = [os.path.abspath(_) for _ in args.reads]
        out_gtf = os.path.abspath(args.out_gtf)
        n_thread = args.threads
        max_indel = args.max_indel
        max_clip = args.max_clip
        min_idt = args.min_idt

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
        input_gtfs = list()
        for i, raw_input_gtf in enumerate(raw_input_gtfs, start=1):
            input_gtf = f"raw_{i}.gtf"
            if os.path.lexists(input_gtf):
                input_gtf = f"raw_{i}_{os.getpid()}.gtf"
            input_gtfs.append(input_gtf)
            os.symlink(raw_input_gtf, input_gtf)

        genome_file = "genome.fasta"
        if os.path.lexists(genome_file):
            genome_file = f"genome_{os.getpid()}.fasta"
        os.symlink(raw_genome_file, genome_file)

        ref_gtf = "ref.gtf"
        if os.path.lexists(ref_gtf):
            ref_gtf = f"ref_{os.getpid()}.gtf"
        os.symlink(raw_ref_gtf, ref_gtf)

        reads_files = list()
        for i, raw_reads_file in enumerate(raw_reads_files, start=1):
            reads_ext = os.path.basename(raw_reads_file).split(".", 1)[1]
            reads_file = f"reads_{i}.{reads_ext}"
            if os.path.lexists(reads_file):
                reads_file = f"reads_{i}_{os.getpid()}.{reads_ext}"
            reads_files.append(reads_file)
            os.symlink(raw_reads_file, reads_file)

        # Main process start
        logger.info(f"Starting isorefiner refine")

        ## Merge step
        run_command(f"gffcompare -o merge -p cons -r {ref_gtf} {ref_gtf} {' '.join(input_gtfs)}", stdout="gffcompare_merge.stdout", stderr="gffcompare_merge.stderr")

        ## Strand-correction step
        run_command(f"gffcompare -r {ref_gtf} merge.combined.gtf -o flip", stdout="gffcompare_flip.stdout", stderr="gffcompare_flip.stderr")
        gtf_flip_strand("merge.combined.gtf", "flip.merge.combined.gtf.tmap", "flipped.gtf")

        ## Correction step based on intron distance 
        run_command(f"gffread -w asm.fa -g {genome_file} flipped.gtf")
        run_command(f"minimap2 -ax map-ont --secondary=no -t {n_thread} asm.fa {' '.join(reads_files)} | samtools view -b -F 2308 -", stdout="raw.bam")
        run_command(f"samtools sort -@ {n_thread} -m 2G raw.bam", stdout="sorted.bam", stderr="samtools_sort.stderr")


        logger.info(f"Finished isorefiner refine")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)


@func_with_log
def gtf_flip_strand(in_gtf, tmap_file, out_gtf):
    target_set = set()
    with open(tmap_file) as fin:
        fin.readline()
        for ln in fin:
            f = ln.rstrip("\n").split("\t")
            if f[2] == "s":
                target_set.add(f[4])

    transcript_re = re.compile(r'transcript_id "([^"]*)"')
    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            m = transcript_re.search(f[8])
            if m and (m.group(1) in target_set):
                if f[6] == "+":
                    f[6] = "-"
                elif f[6] == "-":
                    f[6] = "+"
                print("\t".join(f), file=fout)
            else:
                print(ln, end="", file=fout)
