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
filter_pileup_by_gtf = os.path.join(
    conda_prefix, "opt/popoolation-code/basic-pipeline/filter-pileup-by-gtf.pl"
)

if not os.path.exists(filter_pileup_by_gtf):
    logger.info("Popoolation not installed: checking out code with subversion")
    source = "https://svn.code.sf.net/p/popoolation2/code/trunk"
    popoolation2_code = os.path.join(conda_prefix, "opt/popoolation2-code")
    shell("svn checkout {source} " "{popoolation2_code}")

options = snakemake.params.get("options", "")

shell(
    "perl "
    "{filter_pileup_by_gtf} "
    "{options} "
    "--input {snakemake.input.pileup} "
    "--gtf {snakemake.input.gtf} "
    "--output {snakemake.output.pileup} "
)
