#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

memory = ""
if "mem_mb" in snakemake.resources.keys():
    memory = "-Xmx{}M".format(str(snakemake.resources["mem_mb"]))

shell(
    "picard MarkDuplicates "  # Tool and its subcommand
    "{memory} "  # Automatic Xmx java option
    "{snakemake.params} "  # User defined parmeters
    "--INPUT {snakemake.input} "  # Input file
    "--OUTPUT {snakemake.output.bam} "  # Output bam
    "--METRICS_FILE {snakemake.output.metrics} "  # Output metrics
    "{log}"  # Logging
)
