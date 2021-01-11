"""
Wrapper script to run GENECONV
"""

import sys
import os
import logging
import subprocess


def run_geneconv(cfg_path):
    """
    Runs GENECONV on Windows or Unix-based operating systems
    :param cfg_path: Path to the config file
    :return gc_output: the output of GENECONV
    """

    # Check that the OS is valid
    platform = sys.platform
    try:
        platform.startswith("win") or platform.startswith("linux") or sys.platform == "darwin"
    except OSError:
        print("OSError: {} is not supported".format(sys.platform))

    # Path to GENECONV executables
    script_path = os.path.dirname(os.path.abspath(__file__))

    if sys.platform.startswith("win"):
        bin_path = os.path.join(script_path, 'bin/GENECONV/Windows/geneconv.exe')
    else:
        bin_path = os.path.join(script_path, 'bin/GENECONV/Unix/geneconv')

    if not os.path.isfile(bin_path):
        logging.error("No file exists")

    # Run GENECONV
    if platform.startswith("win"):
        gc_output = subprocess.check_output([bin_path, cfg_path], shell=False, stderr=subprocess.DEVNULL)
    else:
        gc_output = subprocess.check_output([bin_path, cfg_path], shell=False, stderr=subprocess.STDOUT)

    return gc_output
