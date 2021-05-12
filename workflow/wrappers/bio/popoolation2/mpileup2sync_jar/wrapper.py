#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
import sys
import subprocess
import tempfile
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)
mem_mb = snakemake.resources.get("mem_mb", 1000)

conda_prefix = os.getenv("CONDA_PREFIX")
mpileup2sync = os.path.join(conda_prefix, "opt/popoolation2-code/mpileup2sync.jar")

options = snakemake.params.get("options", "")
mpileup = snakemake.input.mpileup
gzip = "" if snakemake.wildcards.gz == "" else "| gzip -vc "

tmp = os.path.basename(tempfile.mkstemp()[1])
fifo = f"{mpileup}{tmp}.fifo"

if os.path.exists(fifo):
    os.unlink(fifo)

shell("mkfifo {fifo} {log}")
shell("zcat {mpileup} > {fifo} &")
args = [
    "java",
    f"-Xmx{mem_mb}M",
    "-jar",
    f"{mpileup2sync}",
    f"{options}",
    "--input",
    f"{fifo}",
    "--output",
    "/dev/stdout",
    "--threads",
    f"{snakemake.threads}",
    f"{log}",
    f"{gzip}",
    "-",
    ">",
    f"{snakemake.output.sync}",
    f"{log}",
]
res = subprocess.run(" ".join(args), shell=True, check=False)

if os.path.exists(fifo):
    os.unlink(fifo)
