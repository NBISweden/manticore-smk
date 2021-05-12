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
fst_sliding = os.path.join(conda_prefix, "opt/popoolation2-code/fst-sliding.pl")

options = snakemake.params.get("options", "")
sync = snakemake.input.sync
outfile = os.path.splitext(snakemake.output.fst)[0]


samples = snakemake.params.samples
config = snakemake.config
sex = snakemake.wildcards.sex
if sex not in config["workflow"]["regions"][snakemake.wildcards.region]["ploidy"]:
    sex = "common"
ploidy = config["workflow"]["regions"][snakemake.wildcards.region]["ploidy"][sex]
pool_size = samples["size"].tolist()
if len(set(pool_size)) == 1:
    pool_size = pool_size[0] * ploidy
else:
    pool_size = ":".join([x * ploidy for x in pool_size])

## NB: fifos won't work because the script reads the input twice to get the max coverage
shell("gzip -fkdv {sync} {log}")
syncplain = os.path.splitext(str(sync))[0]
shell(
    "perl {fst_sliding} {options} "
    "--window-size {snakemake.wildcards.window_size} "
    "--step-size {snakemake.wildcards.step_size} "
    "--pool-size {pool_size} "
    "--input {syncplain} "
    "--output {outfile} {log}"
)
shell("gzip -vf {outfile} {log}")
shell("rm -f {syncplain}")
