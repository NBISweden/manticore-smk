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

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
zlog = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

shell(
    "zcat {snakemake.input.sync} | "
    "perl -pe "
    "'$_=~s/(d+:\d+:\d+:\d+:0:)\d+(?)/${{1}}0/g' "
    "{log} | "
    "gzip -v - {zlog} > {snakemake.output.sync}"
)
