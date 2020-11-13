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
zlog = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

options = snakemake.params.get("options", "")
targets = snakemake.input.targets
gz = snakemake.wildcards.gz
gzip = "" if gz == "" else "| gzip -v "


shell(
    "samtools mpileup "
    "-l {targets} "
    "-B {snakemake.input.bam} {log}"
    "{gzip} > "
    "{snakemake.output.pileup}"
    "{options} {zlog}"
)
