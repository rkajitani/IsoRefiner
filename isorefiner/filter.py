#!/usr/bin/python

import logging
import os
import subprocess
import sys
from isorefiner.common import run_command


def main(args):
    try:
        # Set variables
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

        # Set input files and create symbolic links
        input_gtf = "raw.gtf"
        if os.path.lexists(input_gtf):
            input_gtf = f"raw_{os.getpid()}.gtf"
        genome_fasta = "genome.fasta"
        if os.path.lexists(genome_fasta):
            genome_fasta = f"genome_{os.getpid()}.fasta"
        reads_ext = os.path.basename(os.path.abspath(args.reads_fastq)).split(".", 1)[1]
        reads_fastq = f"reads.{reads_ext}"
        if os.path.lexists(reads_fastq):
            reads_fastq = f"reads_{os.getpid()}.{reads_ext}"
        os.symlink(os.path.abspath(args.input_gtf), input_gtf)
        os.symlink(os.path.abspath(args.genome_fasta), genome_fasta)
        os.symlink(os.path.abspath(args.reads_fastq), reads_fastq)

        # Set logger
        logging.basicConfig(
            filename="log.txt",
            level=logging.INFO,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        logger = logging.getLogger(__name__)

        logger.info(f"Starting isorefiner filter (PID: {os.getpid()})")
        run_command(cmd=f"echo 'gffread -w asm.fa -g genome.fa raw.gtf'", logger=logger)

    except Exception as e:
        print("Error:", e, file=sys.stderr)
        print("Exception class:", type(e).__name__, file=sys.stderr)
        sys.exit(1)


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
