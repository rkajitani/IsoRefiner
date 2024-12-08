#!/usr/bin/env python

import logging
import subprocess
from functools import wraps
import pysam


logger = logging.getLogger(__name__)


def func_with_log(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        func_log_str = f"function {f.__name__}, args={args}, kwargs={kwargs}"
        logger.info(f"Starting {func_log_str}")
        v = f(*args, **kwargs)
        logger.info(f"Finished {func_log_str}")
        return v
    return wrapper


def run_command(cmd, stdout=None, stderr=None):
    try:
        base_cmd = cmd.split(" ")[0]
        if stdout is None:
            stdout = f"{base_cmd}.stdout"
        if stderr is None:
            stderr = f"{base_cmd}.stderr"
        cmd_log_str = f"{cmd}, stdout={stdout}, stderr={stderr}"
        logger.info(f"Starting command: {cmd_log_str}")
        with open(stdout, "w") as out_file, open(stderr, "a") as err_file:
            subprocess.run(
                cmd,
                shell=True,
                stdout=out_file,
                stderr=err_file,
                check=True,
                text=True
            )
        logger.info(f"Finished command: {cmd_log_str}")
    except Exception as e:
        raise e


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
