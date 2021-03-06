#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
import tempfile
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
zlog = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

conda_prefix = os.getenv("CONDA_PREFIX")
script = os.path.join(
    conda_prefix, "opt/popoolation-code/basic-pipeline/subsample-pileup.pl"
)

options = snakemake.params.get("options", "")
tmp = os.path.basename(tempfile.mkstemp()[1])
pileup = snakemake.input.pileup
fifo = f"{pileup}{tmp}.fifo"
if os.path.exists(fifo):
    os.unlink(fifo)

output = os.path.splitext(snakemake.output.pileup)[0]

shell("mkfifo {fifo}")
shell("zcat {pileup} > {fifo} &")
shell("perl " "{script} " "{options} " "--input {fifo} " "--output {output}" " {log}")
shell("gzip -v {output} {zlog}")
if os.path.exists(fifo):
    os.unlink(fifo)
