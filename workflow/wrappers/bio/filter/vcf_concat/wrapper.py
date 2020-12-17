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

log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

tool = snakemake.wildcards.tool
options = snakemake.params.get("options", "")
if options is None:
    options = ""
java_opts = snakemake.params.get("java_opts", "")
ref = snakemake.input.ref

fh = open(snakemake.output.cmd, "w")
fh.write("#!/bin/bash\n")

if tool == "bcftools":
    cmd = (
        f"bcftools concat {options} -O z "
        f"{snakemake.input.vcf} {log} | "
        f"bcftools sort -O z {log} | "
        f"tee {snakemake.output.vcf} | "
        f"bcftools stats -F {ref} - > {snakemake.output.stats} {log} "
    )
else:
    logger.error(f"Tool {tool} not defined for filter {snakemake.wildcards.filtername}")

fh.write(cmd)
fh.close()

shell(cmd)
shell("bcftools index -t {snakemake.output.vcf}")
