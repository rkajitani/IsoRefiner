#!/usr/bin/env python

import logging
import subprocess


logger = logging.getLogger(__name__)


def run_command(cmd, stdout=None, stderr=None):
    try:
        base_cmd = cmd.split(" ")[0]
        if stdout is None:
            stdout = f"{base_cmd}.stdout"
        if stderr is None:
            stderr = f"{base_cmd}.stderr"
        logger.info(f"Running command: {cmd}, stdout: {stdout}, stderr: {stderr}")
        with open(stdout, "w") as out_file, open(stderr, "a") as err_file:
            subprocess.run(
                cmd,
                shell=True,
                stdout=out_file,
                stderr=err_file,
                check=True,
                text=True
            )
        logger.info(f"Finished command: {cmd}")
    except Exception as e:
        raise e
