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

source = snakemake.input.uri
if len(source) == 0:
    source = snakemake.params.uri
dest = snakemake.output.resource_file
scheme = snakemake.params.scheme
cmd_map = {"": "ln -s", "rsync": "rsync -av", "sftp": "cp", "file": "ln -s"}
cmd = cmd_map[scheme]


logger.debug(f"Using command {cmd} {source} {dest} {log}")
shell("{cmd} {source} {dest} {log}")
