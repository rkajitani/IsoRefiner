#!/usr/bin/env python

import subprocess


def run_command(cmd, out_file_name="", err_file_name=""):
    try:
        base_cmd = cmd.split(" ")[0]
        if out_file_name == "":
            out_file_name = f"{base_cmd}.stdout"
        if err_file_name == "":
            err_file_name = f"{base_cmd}.stderr"
        with open(out_file_name, "w") as out_file, open(err_file_name, "a") as err_file:
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
