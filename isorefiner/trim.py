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
        raw_reads_files = [os.path.abspath(_) for _ in args.reads]
        out_prefix = os.path.abspath(args.out_prefix)
        n_thread = args.threads
        tool_option = args.tool_option

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

        n_file = len(raw_reads_files)
        out_files = list()
        if n_file == 1:
            out_files.append(f"{out_prefix}.{reads_ext}")
        else:
            for i, raw_reads_file in enumerate(raw_reads_files, start=1):
                reads_ext = os.path.basename(raw_reads_file).split(".", 1)[1]
                out_files.append(f"{out_prefix}_{i}.{reads_ext}")

        # Main process
        logger.info(f"isorefiner version: {get_version()}")
        logger.info(f"Starting isorefiner {args.command}")
        for i, (reads_file, out_file) in enumerate(zip(reads_files, out_files), start=1):
            run_command(
                " ".join([
                    "porechop_abi",
                    "--ab_initio",
                    "--verbosity 1",
                    f"--threads {n_thread}",
                    f"--input {reads_file}",
                    f"--output {out_file}",
                    f"--temp_dir tmp_{i}",
                    tool_option
                ]),
                stdout=f"porechop_abi_{i}.stdout",
                stderr=f"porechop_abi_{i}.stderr"
            )
            shutil.rmtree(f"tmp_{i}")
        logger.info(f"Finished isorefiner {args.command}")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)
