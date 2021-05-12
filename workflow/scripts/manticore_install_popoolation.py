#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import argparse
from snakemake.deployment import conda
from snakemake import workflow, shell
from snakemake.logging import logger

is_script = False

wrapper_prefix = os.path.join(
    os.path.dirname(os.path.dirname(os.path.realpath((__file__)))), "wrappers"
)
envfiles = {
    "popoolation": [
        os.path.join(
            wrapper_prefix,
            "bio",
            "popoolation",
            "filter_pileup_by_gtf",
            "environment.yaml",
        ),
        os.path.join(
            wrapper_prefix, "bio", "popoolation", "subsample_pileup", "environment.yaml"
        ),
        os.path.join(
            wrapper_prefix, "bio", "popoolation", "variance_sliding", "environment.yaml"
        ),
    ],
    "popoolation2": [
        os.path.join(
            wrapper_prefix, "bio", "popoolation2", "fisher_test", "environment.yaml"
        ),
        os.path.join(
            wrapper_prefix, "bio", "popoolation2", "fst_sliding", "environment.yaml"
        ),
        os.path.join(
            wrapper_prefix,
            "bio",
            "popoolation2",
            "indel_filtering_filter_sync_by_gtf",
            "environment.yaml",
        ),
        os.path.join(
            wrapper_prefix,
            "bio",
            "popoolation2",
            "indel_filtering_identify_indel_regions",
            "environment.yaml",
        ),
        os.path.join(
            wrapper_prefix,
            "bio",
            "popoolation2",
            "mpileup2sync_jar",
            "environment.yaml",
        ),
        os.path.join(
            wrapper_prefix,
            "bio",
            "popoolation2",
            "snp_frequency_diff",
            "environment.yaml",
        ),
        os.path.join(
            wrapper_prefix,
            "bio",
            "popoolation2",
            "subsample_synchronized",
            "environment.yaml",
        ),
    ],
}

try:
    from snakemake.script import Snakemake

    is_script = isinstance(snakemake, Snakemake)
except NameError as e:
    pass

if is_script:
    env_dir = os.path.realpath(snakemake.input.env_dir)
    snakefile = os.path.realpath(snakemake.input.snakefile)
else:
    import snakemake
    from snakemake.deployment import conda

    parser = argparse.ArgumentParser(
        description="Install popoolation and popoolation2 in a conda environment"
    )
    parser.add_argument(
        "--snakefile",
        "-s",
        dest="snakefile",
        type=str,
        help="snakefile",
        default="Snakefile",
    )
    parser.add_argument(
        "--wrapper-prefix",
        "-w",
        dest="wrapper_prefix",
        type=str,
        help="wrapper prefix",
        default=wrapper_prefix,
    )
    parser.add_argument(
        "--env-dir",
        "-d",
        dest="env_dir",
        type=str,
        help="conda prefix",
        default=".snakemake/conda",
    )
    args = parser.parse_args()
    env_dir = os.path.realpath(args.env_dir)
    snakefile = os.path.realpath(args.snakefile)

wf = workflow.Workflow(snakefile)
for k, flist in envfiles.items():
    for fn in flist:
        env = conda.Env(fn, wf, env_dir=env_dir)
        conda_prefix = os.path.join(env_dir, env.hash)
        code = os.path.join(conda_prefix, f"opt/{k}-code")
        logger.info(f"Checking if {code} exists")
        if not os.path.exists(code):
            logger.info(f"installing {k}-code to {code}")
            shell("svn checkout https://svn.code.sf.net/p/{k}/code/trunk " "{code}")
