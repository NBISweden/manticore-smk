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

wc = snakemake.wildcards
ploidy = f"-ploidy {config['workflow']['regions'][wc.region]['ploidy']}"

mode = " -ERC GVCF " if snakemake.wildcards.mode == ".g" else " "

gatk = "GATK"
options = snakemake.params.get("options", "")
java_opts = snakemake.params.get("java_opts", "")

bams = snakemake.input.bam
if isinstance(bams, str):
    bams = [bams]
bams = list(map("-I {}".format, bams))

shell(
    "gatk --java-options '{java_opts}' HaplotypeCaller {options} "
    "-R {snakemake.input.ref} {bams} "
    " {mode} {ploidy} "
    "-O {snakemake.output.vcf} {log}"
)