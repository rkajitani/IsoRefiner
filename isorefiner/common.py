#!/usr/bin/env python

import logging
import subprocess
from functools import wraps


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
