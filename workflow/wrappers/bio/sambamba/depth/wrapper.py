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

inbam = snakemake.input.bam
outbed = snakemake.output.bed

threads = snakemake.threads
extra = snakemake.params.get("extra", "-c 0")
gzip = ""
if re.search("(.gz|.gzip)$", snakemake.output.bed):
    gzip = " | gzip -v - "

window_size = snakemake.params.window_size

shell(
    "sambamba depth window -t {threads} {extra} -w {window_size} {inbam}"
    " {gzip} > {outbed} {log}"
)
