#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
import gzip
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

seq = snakemake.input.seq
out = snakemake.output.jf
kmer = snakemake.wildcards.kmer
ext = snakemake.wildcards.ext
options = snakemake.params.options

cat = "cat"

if ext == ".txt":
    if re.search(r".gz$", seq):
        with gzopen(seq, "rb") as fh:
            seqin = fh.readlines()
    else:
        with open(seq, "r") as fh:
            seqin = fh.readlines()
else:
    seqin = seq

seqin = seq

nzip = sum([re.search(r".gz$", x) is None for x in seqin])
assert nzip == 0 or len(seqin) == nzip, "infiles must all be unzipped or gzipped!"

if nzip == len(seqin):
    cat = "zcat"

shell(
    "zcat {seqin} | jellyfish count {options} -m {kmer} -t {snakemake.threads} /dev/fd/0 -o {out} {log}"
)
