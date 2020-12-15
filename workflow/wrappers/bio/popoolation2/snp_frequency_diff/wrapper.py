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
snp_frequency_diff = os.path.join(
    conda_prefix, "opt/popoolation2-code/snp-frequency-diff.pl"
)

if not os.path.exists(snp_frequency_diff):
    logger.info("Popoolation2 not installed: checking out code with subversion")
    popoolation2_code = os.path.join(conda_prefix, "opt/popoolation2-code")
    shell(
        "svn checkout https://svn.code.sf.net/p/popoolation2/code/trunk "
        "{popoolation2_code}"
    )

options = snakemake.params.get("options", "")
sync = snakemake.input.sync
rc_out = re.sub(".gz$", "", snakemake.output.rc)
pwc_out = re.sub(".gz$", "", snakemake.output.pwc)
outprefix = re.sub("_rc.txt$", "", rc_out)

shell("gzip -fkdv {sync} {log}")
syncplain = re.sub(".gz$", "", sync)
shell(
    "perl {snp_frequency_diff} {options} "
    "--input {syncplain} "
    "--output-prefix {outprefix} {log}"
)
shell("gzip -vcf {rc_out} {log} > {snakemake.output.rc}")
shell("gzip -vcf {pwc_out} {log} > {snakemake.output.pwc}")
shell("rm -f {syncplain}")
