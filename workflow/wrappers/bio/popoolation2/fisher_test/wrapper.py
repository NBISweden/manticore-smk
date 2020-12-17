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

conda_prefix = os.getenv("CONDA_PREFIX")
fisher_test = os.path.join(conda_prefix, "opt/popoolation2-code/fisher-test.pl")

if not os.path.exists(fisher_test):
    logger.info("Popoolation2 not installed: checking out code with subversion")
    source = "https://svn.code.sf.net/p/popoolation2/code/trunk"
    popoolation2_code = os.path.join(conda_prefix, "opt/popoolation2-code")
    shell("svn checkout {source} " "{popoolation2_code}")

# Add Text/NSP/Measures/2D/Fisher/twotailed.pm to search path
siteperl = f"{conda_prefix}/lib/perl5/site_perl"
for v in os.listdir(siteperl):
    os.environ["PERL5LIB"] = f"{siteperl}/{v}:{os.environ.get('PERL5LIB', '')}"

options = snakemake.params.get("options", "")
sync = snakemake.input.sync
outfile = os.path.splitext(snakemake.output.fet)[0]

shell("gzip -fkdv {sync} {log}")
syncplain = os.path.splitext(str(sync))[0]
shell(
    "perl {fisher_test} {options} "
    "--window-size {snakemake.wildcards.window_size} "
    "--step-size {snakemake.wildcards.step_size} "
    "--input {syncplain} "
    "--output {outfile} {log}"
)
shell("gzip -vf {outfile} {log}")
shell("rm -f {syncplain}")
