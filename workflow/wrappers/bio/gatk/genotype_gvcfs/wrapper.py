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
annotation = snakemake.params.get("annotation", [])

ref = f"-R {snakemake.input.ref} "
db = "-V gendb://{}".format(snakemake.input.db)
annotation_opt = [f"--annotation {a}" for a in annotation]


shell(
    "gatk --java-options '{java_opts}' GenotypeGVCFs {options} "
    "{ref} {targets} {annotation_opt} --include-non-variant-sites true "
    "-O {snakemake.output.vcf} {db} {log}"
)
