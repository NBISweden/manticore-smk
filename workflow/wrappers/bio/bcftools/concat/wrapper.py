#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)
options = snakemake.params.get("extra", "")
ref = snakemake.input.ref

shell(
    "bcftools concat {options} -O z "
    "{snakemake.input.vcf} {log} | "
    "bcftools sort -O z {log} | "
    "tee {snakemake.output.vcf} | "
    "bcftools stats -F {ref} - > {snakemake.output.stats}"
    " {log}"
)
shell("bcftools index -t {snakemake.output.vcf}")
