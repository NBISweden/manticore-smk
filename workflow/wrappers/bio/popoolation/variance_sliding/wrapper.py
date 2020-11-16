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

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
zlog = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

conda_prefix = os.getenv("CONDA_PREFIX")
script = os.path.join(conda_prefix, "opt/popoolation-code/Variance-sliding.pl")

if not os.path.exists(script):
    logger.info("Popoolation not installed: checking out code with subversion")
    popoolation_code = os.path.join(conda_prefix, "opt/popoolation-code")
    shell(
        "svn checkout https://svn.code.sf.net/p/popoolation/code/trunk "
        "{popoolation_code}"
    )

options = snakemake.params.get("options", "")
samples = snakemake.params.samples
config = snakemake.config
sex = snakemake.wildcards.get("sex", "common")
pool_size = samples[samples.SM == snakemake.wildcards.sample]
pool_size = (
    pool_size["size"].to_list()[0]
    * config["workflow"]["regions"][snakemake.wildcards.region]["ploidy"][sex]
)
outtxt = os.path.splitext(snakemake.output.txt)[0]

shell(
    "perl "
    "{script} "
    "{options} "
    "--measure {snakemake.wildcards.measure} "
    "--pool-size {pool_size} "
    "--input {snakemake.input.pileup} "
    "--output {outtxt} "
    "--window-size {snakemake.wildcards.window_size} "
    "--step-size {snakemake.wildcards.step_size} "
    "{log} && gzip -v {outtxt} {zlog}"
)
