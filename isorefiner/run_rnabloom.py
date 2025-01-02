#!/usr/bin/python

import logging
import os
import sys
from isorefiner.common import run_command, func_with_log

logger = logging.getLogger(__name__)

        
def main(args):
    try:
        # Set variables
        raw_reads_files = [os.path.abspath(_) for _ in args.reads]
        raw_genome_file = os.path.abspath(args.genome)
        out_gtf = os.path.abspath(args.out_gtf)
        n_thread = args.threads
        max_mem = args.max_mem
        rnabloom_option = args.rnabloom_option
        gmap_min_cov = args.gmap_min_cov
        gmap_min_idt = args.gmap_min_idt
        gmap_max_intron = args.gmap_max_intron
        args.gmap_option = args.gmap_option

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
        reads_files = list()
        for i, raw_reads_file in enumerate(raw_reads_files, start=1):
            reads_ext = os.path.basename(raw_reads_file).split(".", 1)[1]
            reads_file = f"reads_{i}.{reads_ext}"
            if os.path.lexists(reads_file):
                reads_file = f"reads_{i}_{os.getpid()}.{reads_ext}"
            reads_files.append(reads_file)
            os.symlink(raw_reads_file, reads_file)

        # Main process
        logger.info(f"Starting isorefiner {args.command}")
        ## RNA-Bloom
        os.environ["JAVA_TOOL_OPTIONS"] = f"-Xmx{max_mem}"
        run_command(f"rnabloom -t {n_thread} -outdir rnabloom_out {rnabloom_option} -long {' '.join(reads_files)}")
        ## GMAP
        logger.info(f"Finished isorefiner {args.command}")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)
