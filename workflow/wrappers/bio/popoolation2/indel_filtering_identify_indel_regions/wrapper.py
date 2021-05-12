#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
import tempfile
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

conda_prefix = os.getenv("CONDA_PREFIX")
script = os.path.join(
    conda_prefix, "opt/popoolation2-code/indel_filtering/identify-indel-regions.pl"
)

options = snakemake.params.get("options", "")

mpileup = snakemake.input.mpileup
tmp = os.path.basename(tempfile.mkstemp()[1])
fifo = f"{mpileup}{tmp}.fifo"
if os.path.exists(fifo):
    os.unlink(fifo)

shell("mkfifo {fifo}")
shell("zcat {mpileup} > {fifo} &")
shell(
    "perl "
    "{script} "
    "{options} "
    "--input {fifo} "
    "--output {snakemake.output.gtf} "
    "{log}"
)
if os.path.exists(fifo):
    os.unlink(fifo)
