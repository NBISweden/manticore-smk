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
config = snakemake.config

targets = snakemake.input.get("targets", "")
targets = f" -L {targets}" if targets else ""


gatk = "GATK"
options = snakemake.params.get("options", "")
java_opts = snakemake.params.get("java_opts", "")
output = snakemake.output

vcfs = snakemake.input.vcf
vcfs = list(map("-V {}".format, vcfs))

if "mem_mb" in snakemake.resources.keys():
    java_opts += " -Xmx{}M".format(str(snakemake.resources["mem_mb"]))

shell(
    "gatk --java-options '{java_opts}' GenomicsDBImport {options} "
    "-R {snakemake.input.ref} {targets} {vcfs} "
    "--genomicsdb-workspace-path {output} "
    " {log} "
)
