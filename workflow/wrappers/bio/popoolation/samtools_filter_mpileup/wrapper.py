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

options = snakemake.params.get("options", {})
filter_options = options.get("filter_options", "")
mpileup_options = options.get("mpileup_options", "")

targets = snakemake.input.targets
gz = snakemake.wildcards.gz
gzip = "" if gz == "" else "| gzip -v "

shell(
    "samtools view "
    "-L {targets} "
    "{filter_options} "
    "{snakemake.input.bam} "
    " -b | samtools mpileup {mpileup_options} - "
    "{gzip} "
    "> {snakemake.output.pileup}"
)
