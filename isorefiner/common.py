#!/usr/bin/env python

import subprocess


def run_command(cmd, out_file_name=None, err_file_name=None, logger=None):
    try:
        base_cmd = cmd.split(" ")[0]
        if out_file_name is None:
            out_file_name = f"{base_cmd}.stdout"
        if err_file_name is None:
            err_file_name = f"{base_cmd}.stderr"
        if logger is not None:
            logger.info(f"Running command: {cmd}, stdout: {out_file_name}, stderr: {err_file_name}")
        with open(out_file_name, "w") as out_file, open(err_file_name, "a") as err_file:
            subprocess.run(
                cmd,
                shell=True,
                stdout=out_file,
                stderr=err_file,
                check=True,
                text=True
            )
        if logger is not None:
            logger.info(f"Finished command: {cmd}")
    except Exception as e:
        raise e
