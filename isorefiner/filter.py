#!/usr/bin/python

import logging
import os
import subprocess
import sys
import polars as pl
import pysam
from isorefiner.common import func_with_log, run_command

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
        min_cov = args.min_cov
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
        input_gtf = "raw.gtf"
        if os.path.lexists(input_gtf):
            input_gtf = f"raw_{os.getpid()}.gtf"
        os.symlink(raw_input_gtf, input_gtf)

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
        run_command(f"gffread -w asm.fa -g {genome_file} {input_gtf}")
        run_command(f"minimap2 -ax map-ont --secondary=no -t {n_thread} asm.fa {' '.join(reads_files)} | samtools view -b -F 2308 -", stdout="raw.bam")
        run_command(f"samtools sort -@ {n_thread} -m 2G raw.bam", stdout="sorted.bam", stderr="samtools_sort.stderr")
        filter_bam("sorted.bam", "filt.bam", max_indel, max_clip, min_idt)
        run_command(f"samtools index filt.bam")
        run_command(f"samtools coverage filt.bam", stdout="cov.tsv")
        filter_coverage_table("cov.tsv", "cov_depth_filt.list", min_cov, min_mean_depth)

        logger.info(f"Finished isorefiner filter")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)


@func_with_log
def filter_bam(in_bam, out_bam, max_indel, max_clip, min_idt):
    with pysam.AlignmentFile(in_bam, "rb") as fin, pysam.AlignmentFile(out_bam, "wb", header=fin.header) as fout:
        for aln in fin:
            valid_flag = True
            idt = 1.0 - (aln.get_tag("NM") / aln.query_alignment_length)
            if idt < min_idt:
                continue
            for cig_ope, cig_len in aln.cigartuples:
                if (((cig_ope == 1 or cig_ope == 2) and cig_len > max_indel) or
                    (cig_ope == 4 or cig_ope == 5) and cig_len > max_clip):
                    valid_flag = False
                    break
            if valid_flag:
                fout.write(aln)


@func_with_log
def filter_coverage_table(in_file, out_file, min_cov, min_mean_depth):
    df = (pl.read_csv(in_file, separator="\t")
        .select("#rname", "coverage", "meandepth")
        .filter((pl.col("coverage") >= min_cov) & (pl.col("meandepth") >= min_mean_depth))
        .select("#rname")
    )
    df.write_csv(out_file, include_header=False)


#../test/bambu/filter/cmd.sh
#
#t=4
#max_indel=20
#max_clip=200
#min_idt=0.9
#min_cov=95
#min_mean_depth=1
#
#ln -s ../out_transcript.gtf raw.gtf
#ln -s ../../reads.fastq.gz
#ln -s ../../genome_chr22.fa genome.fa
#
#gffread -w asm.fa -g genome.fa raw.gtf
#/usr/bin/time minimap2 -ax map-ont --secondary=no -t $t asm.fa reads.fastq.gz 2>mm2.log | samtools view -b -F 2308 - >raw.bam
#samtools sort -@ $t -m 2G raw.bam >sorted.bam 2>sort.log
#./bam_filter.py sorted.bam filt.bam $max_indel $max_clip $min_idt
#samtools index filt.bam
#samtools coverage filt.bam >cov.tsv
#
#df_select.py cov.tsv "#rname" coverage meandepth | sed 1d | perl -ane "print(\$F[0], \"\n\") if (\$F[1] >= $min_cov and \$F[2] >= $min_mean_depth)" >cov_depth_filt.list
#./gtf_extract_transcript.py raw.gtf cov_depth_filt.list >cov_depth_filt.gtf
#ln -s cov_depth_filt.gtf out_transcript.gtf
