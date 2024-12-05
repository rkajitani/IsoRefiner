#!/usr/bin/env python

import os
import subprocess


def run_command(cmd, out_file_name, err_file_name):
    try:
        with open(out_file_name, "w") as out_file, open(err_file_name, "w") as err_file:
            subprocess.run(
                cmd,
                shell=True,
                stdout=out_file,
                stderr=err_file,
                check=True,
                text=True
            )
    except Exception as e:
        raise e
