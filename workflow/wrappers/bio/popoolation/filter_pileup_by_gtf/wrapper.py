#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

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
pileup = snakemake.output.pileup

# FIXME: output is truncated if zipped input/output
shell(
    "perl "
    "{script} "
    "{options} "
    "--input {snakemake.input.pileup} "
    "--gtf {snakemake.input.gtf} "
    "--output {pileup} {log}"
)
