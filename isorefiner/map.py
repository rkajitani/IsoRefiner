#!/usr/bin/python

import logging
import os
import sys
from isorefiner.common import func_with_log, run_command

logger = logging.getLogger(__name__)

        
def main(args):
    try:
        # Set variables
        raw_reads_files = [os.path.abspath(_) for _ in args.reads]
        raw_genome_file = os.path.abspath(args.genome)
        out_prefix = os.path.abspath(args.out_prefix)
        n_thread = args.threads
        mm2_option = args.mm2_option
        sort_option = args.sort_option

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
        reads_files = list()
        for i, raw_reads_file in enumerate(raw_reads_files, start=1):
            reads_ext = os.path.basename(raw_reads_file).split(".", 1)[1]
            reads_file = f"reads_{i}.{reads_ext}"
            if os.path.lexists(reads_file):
                reads_file = f"reads_{i}_{os.getpid()}.{reads_ext}"
            reads_files.append(reads_file)
            os.symlink(raw_reads_file, reads_file)

        n_file = len(raw_reads_files)
        out_files = list()
        if n_file == 1:
            out_files.append(f"{out_prefix}.bam")
        else:
            for i, raw_reads_file in enumerate(raw_reads_files, start=1):
                out_file = f"{out_prefix}_{i}.bam"
                out_files.append(reads_file)

        # Main process
        logger.info(f"Starting isorefiner {args.command}")
        for i, (reads_file, out_file) in enumerate(zip(reads_files, out_files), start=1):
            run_command(
                " ".join([
                    "minimap2",
                    "-a",
                    f"-t {n_thread}",
                    mm2_option,
                    genome_file,
                    reads_file,
                    "| samtools view -Sb"
                ]),
                stdout=f"raw_{i}.bam",
                stderr=f"minimap2_{i}.stderr"
            )
        for i, out_file in enumerate(out_files, start=1):
            run_command(
                " ".join([
                    "samtools sort",
                    "-O BAM",
                    f"-@ {n_thread}",
                    sort_option,
                    ' '.join(reads_files)
                ]),
                stdout=out_file,
                stderr=f"samtools_sort_{i}.stderr"
            )
        for i, out_file in enumerate(out_files, start=1):
            run_command(
                f"samtools index {out_file}",
                stdout=f"samtools_index_{i}.stdout",
                stderr=f"samtools_index_{i}.stderr"
            )
        logger.info(f"Finished isorefiner {args.command}")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)
