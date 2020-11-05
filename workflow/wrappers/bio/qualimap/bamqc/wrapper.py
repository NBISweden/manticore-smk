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

bam = snakemake.input.bam
outdir = os.path.dirname(snakemake.output.genome_results)
mem_mb = snakemake.resources.mem_mb
threads = snakemake.threads

shell(
    "unset DISPLAY; "
    "qualimap bamqc -bam {bam} -nt {threads} qualimap --java-mem-size={mem_mb}M -outdir {outdir}"
    " {log}"
)
