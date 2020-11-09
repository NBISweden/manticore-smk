#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

from snakemake.shell import shell


inputs = " ".join(f"--INPUT {f}" for f in snakemake.input.vcfs)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

shell(
    "picard"
    " MergeVcfs"
    " {extra}"
    " {inputs}"
    " --OUTPUT {snakemake.output[0]}"
    " {log}"
)
