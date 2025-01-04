#!/usr/bin/python

import logging
import os
import shutil
import sys
from isorefiner.common import get_version, run_command

logger = logging.getLogger(__name__)

        
def main(args):
    try:
        # Set variables
        raw_bam_files = [os.path.abspath(_) for _ in args.bam]
        raw_genome_file = os.path.abspath(args.genome)
        raw_ref_gtf = os.path.abspath(args.ref_gtf)
        out_gtf = os.path.abspath(args.out_gtf)
        n_thread = args.threads
        tool_s_option = args.tool_s_option
        tool_c_option = args.tool_c_option
        tool_q_option = args.tool_q_option

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
        with open("samples.tsv", "w") as fout:
            for i, bam_file in enumerate(bam_files):
                print(f"{bam_file}\tsample_{i}", file=fout)

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
        run_command(
            " ".join([
                "perl $(which ESPRESSO_S.pl)",
                f"-T {n_thread}",
                "-L samples.tsv",
                f"-F {genome_file}",
                f"-A {ref_gtf}",
                "-O espresso_work",
                tool_s_option
            ]),
            stdout=f"ESPRESSO_S.stdout",
            stderr=f"ESPRESSO_S.stderr"
        )
        for i in range(0, len(bam_files)):
            run_command(
                " ".join([
                    "perl $(which ESPRESSO_C.pl)",
                    f"-X {i}",
                    f"-T {max(n_thread, 2)}",
                    "-I espresso_work",
                    f"-F {genome_file}",
                    tool_c_option
                ]),
                stdout=f"ESPRESSO_C_{i}.stdout",
                stderr=f"ESPRESSO_C_{i}.stderr"
            )
        run_command(
            " ".join([
                "perl $(which ESPRESSO_Q.pl)",
                f"-T {n_thread}",
                "-L espresso_work/samples.tsv.updated",
                f"-A {ref_gtf}",
                tool_q_option
            ]),
            stdout=f"ESPRESSO_Q.stdout",
            stderr=f"ESPRESSO_Q.stderr"
        )
        shutil.move("espresso_work/samples_N2_R0_updated.gtf", out_gtf)
        logger.info(f"Finished isorefiner {args.command}")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)
