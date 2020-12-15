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

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

conda_prefix = os.getenv("CONDA_PREFIX")
script = os.path.join(
    conda_prefix, "opt/popoolation-code/basic-pipeline/filter-pileup-by-gtf.pl"
)

if not os.path.exists(script):
    logger.info("Popoolation not installed: checking out code with subversion")
    popoolation_code = os.path.join(conda_prefix, "opt/popoolation-code")
    shell(
        "svn checkout https://svn.code.sf.net/p/popoolation/code/trunk "
        "{popoolation_code}"
    )

options = snakemake.params.get("options", "")
pileup = snakemake.input.pileup
outfile = re.sub(".gz$", "", str(snakemake.output.pileup))
tmp = os.path.basename(tempfile.mkstemp()[1])
fifo = f"{pileup}{tmp}.fifo"
if os.path.exists(fifo):
    os.unlink(fifo)

shell("mkfifo {fifo}")
shell("zcat {pileup} > {fifo} &")

shell(
    "perl "
    "{script} "
    "{options} "
    "--input {fifo} "
    "--gtf {snakemake.input.gtf} "
    "--output {outfile} {log}"
)
shell("gzip -v {outfile} {log}")

if os.path.exists(fifo):
    os.unlink(fifo)
