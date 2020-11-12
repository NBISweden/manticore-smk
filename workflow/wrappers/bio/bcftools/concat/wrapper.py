#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
options = snakemake.params.get("extra", "")

shell(
    "bcftools concat {options} -O z "
    "{snakemake.input.vcfs} | "
    "tee {snakemake.output.vcf} | "
    "bcftools stats - > {snakemake.output.stats}"
    " {log}"
)
shell("bcftools index -t {snakemake.output.vcf}")
