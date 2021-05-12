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

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

conda_prefix = os.getenv("CONDA_PREFIX")
subsample_sync = os.path.join(
    conda_prefix, "opt/popoolation2-code/subsample-synchronized.pl"
)

options = snakemake.params.get("options", "")
sync = snakemake.input.sync
outfile = os.path.splitext(snakemake.output.sync)[0]

## NB: fifos won't work because the script reads the input twice to get the max coverage
shell("gzip -fkdv {sync} {log}")
syncplain = os.path.splitext(str(sync))[0]

shell("perl {subsample_sync} {options} --input {syncplain} --output {outfile} {log}")
shell("gzip -v {outfile} {log}")
shell("rm -f {syncplain}")
