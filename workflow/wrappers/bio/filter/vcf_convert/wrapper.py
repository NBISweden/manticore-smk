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

fh = open(snakemake.output.cmd, "w")
fh.write("#!/bin/bash\n")

if tool == "gatk":
    cmd = (
        f"gatk --java-options '{java_opts}' SelectVariants -R {snakemake.input.ref} "
        f"-V {snakemake.input.vcf} "
        f"{options} -O {snakemake.output.vcf} {log}"
    )
elif tool == "bcftools":
    cmd = (
        f"bcftools view {options} {snakemake.input.vcf} -O z -o {snakemake.output.vcf} {log};  "
        f"bcftools index --tbi {snakemake.output.vcf} {log}"
    )
else:
    logger.error(f"Tool {tool} not defined for filter {snakemake.wildcards.filtername}")

fh.write(cmd)
fh.close()

shell(cmd)
