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
sex = samples.at[snakemake.wildcards.sample, "sex"]
size = samples.at[snakemake.wildcards.sample, "size"]
if sex not in config["workflow"]["regions"][snakemake.wildcards.region]["ploidy"]:
    sex = "common"
ploidy = config["workflow"]["regions"][snakemake.wildcards.region]["ploidy"][sex]
pool_size = size * ploidy

outtxt = os.path.splitext(snakemake.output.txt)[0]

tmp = os.path.basename(tempfile.mkstemp()[1])
pileup = snakemake.input.pileup
fifo = f"{pileup}{tmp}.fifo"
if os.path.exists(fifo):
    os.unlink(fifo)

shell("mkfifo {fifo}")
shell("zcat {pileup} > {fifo} &")
shell(
    "perl "
    "{script} "
    "{options} "
    "--measure {snakemake.wildcards.measure} "
    "--pool-size {pool_size} "
    "--input {fifo} "
    "--output {outtxt} "
    "--window-size {snakemake.wildcards.window_size} "
    "--step-size {snakemake.wildcards.step_size} "
    "{log} && gzip -v {outtxt} {zlog}"
)
if os.path.exists(fifo):
    os.unlink(fifo)
