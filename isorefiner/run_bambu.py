#!/usr/bin/python

import logging
import os
import shutil
import sys
from isorefiner.common import get_version, run_command, func_with_log

logger = logging.getLogger(__name__)

        
def main(args):
    try:
        # Set variables
        raw_bam_files = [os.path.abspath(_) for _ in args.bam]
        raw_genome_file = os.path.abspath(args.genome)
        raw_ref_gtf = os.path.abspath(args.ref_gtf)
        out_gtf = os.path.abspath(args.out_gtf)
        n_thread = args.threads

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
        bam_files = list()
        for i, raw_bam_file in enumerate(raw_bam_files, start=1):
            bam_file = f"map_{i}.bam"
            if os.path.lexists(bam_file):
                bam_file = f"map_{i}_{os.getpid()}.bam"
            bam_files.append(bam_file)
            os.symlink(raw_bam_file, bam_file)
            for idx_ext in ["bai", "csi"]:
                if os.path.lexists(f"{raw_bam_file}.{idx_ext}"):
                    os.symlink(f"{raw_bam_file}.{idx_ext}", f"{bam_file}.{idx_ext}")
        with open("input_bam.txt", "w") as fout:
            print("\n".join(["#isoquant"] + bam_files), file=fout)

        genome_file = "genome.fasta"
        if os.path.lexists(genome_file):
            genome_file = f"genome_{os.getpid()}.fasta"
        os.symlink(raw_genome_file, genome_file)

        ref_gtf = "ref.gtf"
        if os.path.lexists(ref_gtf):
            ref_gtf = f"ref_{os.getpid()}.gtf"
        os.symlink(raw_ref_gtf, ref_gtf)

        # Main process
        logger.info(f"isorefiner version: {get_version()}")
        logger.info(f"Starting isorefiner {args.command}")
        for bam_file in bam_files:
            if not (os.path.lexists(f"{bam_file}.bai") or os.path.lexists(f"{bam_file}.csi")):
                run_command(f"samtools index {bam_file}.{idx_ext}")
        make_bambu_script("bambu.R")
        run_command(
            f"Rscript --slave --vanilla bambu.R {n_thread} {genome_file} {ref_gtf} {' '.join(bam_files)}",
            stdout=f"bambu.stdout",
            stderr=f"bambu.stderr"
        )
        gtf_filt_zero_count("bambu_out/extended_annotations.gtf", out_gtf, "bambu_out/counts_transcript.txt")
        logger.info(f"Finished isorefiner {args.command}")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)


@func_with_log
def make_bambu_script(script_file):
    with open(script_file, "w") as fout:
        print(r"""
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("usage: Rscript bambu.R n_threads genome.fa annot.gtf bam1 [bam2 ...]\n")
	quit("no", status = 1, runLast = FALSE)
}

library(bambu)

se <- bambu(reads = args[4:length(args)], ncore = args[1], genome = args[2], annotations = args[3])
writeBambuOutput(se, path = "./bambu_out/")
        """.strip(), file=fout)


@func_with_log
def gtf_filt_zero_count(in_gtf, out_gtf, bambu_count_tsv):
    def parse_gtf_attr(attr_field):
        gene_id = ""
        transcript_id = ""
        for attr_str in attr_field.split(";"):
            attr = attr_str.strip().split(" ")
            if len(attr) != 2:
                continue
            if attr[0] == "gene_id":
                gene_id = attr[1].strip('"')
            elif attr[0] == "transcript_id":
                transcript_id = attr[1].strip('"')
        return gene_id, transcript_id

    positive_set = set()
    with open(bambu_count_tsv) as fin:
        fin.readline()
        for ln in fin:
            f = ln.rstrip("\n").split("\t")
            if float(f[2]) > 0.0:
                positive_set.add(f[0])

    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        fin.readline()
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                print(ln, end="", file=fout)
                continue
            f = ln.rstrip("\n").split("\t")
            gene_id, transcript_id = parse_gtf_attr(f[8])
            if (gene_id != "" and
                transcript_id != "" and
                transcript_id in positive_set):
                print(ln, end="", file=fout)