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

options = snakemake.params.get("options", {})
filter_options = options.get("filter_options", "")
mpileup_options = options.get("mpileup_options", "")

targets = snakemake.input.targets
gzip = "" if snakemake.wildcards.gz == "" else "| gzip -v "

shell(
    "samtools view "
    "-L {targets} "
    "{filter_options} "
    "{snakemake.input.bam} "
    " -b | samtools mpileup -a {mpileup_options} - {log}"
    "{gzip} "
    "> {snakemake.output.pileup} {zlog}"
)
