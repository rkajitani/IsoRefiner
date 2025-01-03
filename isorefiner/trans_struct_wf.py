#!/usr/bin/python

import logging
import os
import sys
from isorefiner.common import run_command

logger = logging.getLogger(__name__)

        
def main(args):
    try:
        # Set constant variables
        prefix = "isorefiner"
        mapping_based_tools = ["bambu", "espresso", "isoquant", "stringtie"]
        denovo_tools = ["rnabloom"]

        # Set variables from arguments
        raw_reads_files = [os.path.abspath(_) for _ in args.reads]
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
        genome_file = "genome.fasta"
        if os.path.lexists(genome_file):
            genome_file = f"genome_{os.getpid()}.fasta"
        os.symlink(raw_genome_file, genome_file)

        ref_gtf = "ref.gtf"
        if os.path.lexists(ref_gtf):
            ref_gtf = f"ref_{os.getpid()}.gtf"
        os.symlink(raw_ref_gtf, ref_gtf)

        reads_files = list()
        trimmed_files = list()
        bam_files = list()
        for i, raw_reads_file in enumerate(raw_reads_files, start=1):
            reads_ext = os.path.basename(raw_reads_file).split(".", 1)[1]
            reads_file = f"reads_{i}.{reads_ext}"
            if os.path.lexists(reads_file):
                reads_file = f"reads_{i}_{os.getpid()}.{reads_ext}"
            reads_files.append(reads_file)
            os.symlink(raw_reads_file, reads_file)
            if len(raw_reads_files) > 1:
                trimmed_files.append(f"{prefix}_trimmed_{i}.{reads_ext}")
                bam_files.append(f"{prefix}_mapped_{i}.bam")
            else:
                trimmed_files.append(f"{prefix}_trimmed.{reads_ext}")
                bam_files.append(f"{prefix}_mapped.bam")
        reads_files_str = " ".join(reads_files)
        trimmed_files_str = " ".join(trimmed_files)
        bam_files_str = " ".join(bam_files)

        # Main process
        logger.info(f"Starting isorefiner {args.command}")
        run_command(
            f"isorefiner trim -t {n_thread} -r {reads_files_str}",
            stdout="trim.stdout",
            stderr="trim.stderr"
        )
        run_command(
            f"isorefiner map -t {n_thread} -r {trimmed_files_str} -g {genome_file}",
            stdout="map.stdout",
            stderr="map.stderr"
        )
        for tool in mapping_based_tools:
            run_command(
                f"isorefiner run_{tool} -t {n_thread} -b {bam_files_str} -g {genome_file} -a {ref_gtf}",
                stdout=f"{tool}.stdout",
                stderr=f"{tool}.stderr"
            )
            run_command(
                f"isorefiner filter -t {n_thread} -i {prefix}_{tool}.gtf -r {trimmed_files_str} -g {genome_file} -o {prefix}_filtered_{tool}.gtf -d {prefix}_filter_{tool}_work",
                stdout=f"filter_{tool}.stdout",
                stderr=f"filter_{tool}.stderr"
            )
        for tool in denovo_tools:
            run_command(
                f"isorefiner run_{tool} -t {n_thread} -r {trimmed_files_str} -g {genome_file}",
                stdout=f"{tool}.stdout",
                stderr=f"{tool}.stderr"
            )
            run_command(
                f"isorefiner filter -t {n_thread} -i {prefix}_{tool}.gtf -r {trimmed_files_str} -g {genome_file} -o {prefix}_filtered_{tool}.gtf -d {prefix}_filter_{tool}_work",
                stdout=f"filter_{tool}.stdout",
                stderr=f"filter_{tool}.stderr"
            )
        run_command(
            f"isorefiner refine -t {n_thread} -i {' '.join([f'{prefix}_filtered_{_}.gtf' for _ in mapping_based_tools + denovo_tools])} -r {trimmed_files_str} -g {genome_file} -a {ref_gtf} -o {out_gtf}",
            stdout="refine.stdout",
            stderr="refine.stderr"
        )
        logger.info(f"Finished isorefiner {args.command}")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)
