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
filter_pileup_by_gtf = os.path.join(
    conda_prefix, "opt/popoolation2-code/indel_filtering/filter-sync-by-gtf.pl"
)

if not os.path.exists(filter_pileup_by_gtf):
    logger.info("Popoolation not installed: checking out code with subversion")
    popoolation2_code = os.path.join(conda_prefix, "opt/popoolation2-code")
    shell(
        "svn checkout https://svn.code.sf.net/p/popoolation2/code/trunk "
        "{popoolation2_code}"
    )

options = snakemake.params.get("options", "")
sync = snakemake.input.sync
outfile = os.path.splitext(snakemake.output.sync)[0]

tmp = os.path.basename(tempfile.mkstemp()[1])
fifo = f"{sync}{tmp}.fifo"
if os.path.exists(fifo):
    os.unlink(fifo)

shell("mkfifo {fifo}")
shell("zcat {sync} > {fifo} &")

shell(
    "perl "
    "{filter_pileup_by_gtf} "
    "{options} "
    "--input {fifo} "
    "--gtf {snakemake.input.gtf} "
    "--output {outfile} "
    "{log}"
)
shell("gzip -v {outfile} {log}")

if os.path.exists(fifo):
    os.unlink(fifo)
