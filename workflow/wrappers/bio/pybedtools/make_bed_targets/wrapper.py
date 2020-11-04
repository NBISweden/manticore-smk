#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
import pybedtools
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
regions = list(pybedtools.BedTool(snakemake.input.bed))

p = int(snakemake.wildcards.partition)
npart = int(snakemake.params.npart)

try:
    assert len(regions) >= npart
except AssertionError as e:
    logger.warning(
        (
            f"Number of regions smaller than number of partitions: '{len(regions)} < {npart}': "
            f"lower the number of partitions "
            f"(config['workflow']['regions']['{snakemake.wildcards.region}']['npart'])"
        )
    )
    raise
try:
    assert npart >= p
except AssertionError as e:
    logger.warning(
        (f"Requested partition is larger than number of partitions: '{p} > {npart}'")
    )
    raise

ix = sorted(range(len(regions)), key=lambda k: len(regions[k]), reverse=True)
out = [[regions[i]] for i in ix[0:npart]]
outlen = [len(r[0]) for r in out]
for j in ix[npart : len(ix)]:
    imin = outlen.index(min(outlen))
    out[imin].append(regions[j])
    outlen[imin] += len(regions[j])
bedout = pybedtools.BedTool(
    "\n".join(str(r) for r in out[int(p) - 1]), from_string=True
)
bedout.saveas(snakemake.output.bed)
