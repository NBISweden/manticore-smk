import sys
import re
import os
import itertools
import urllib
import copy
import contextlib
import subprocess as sp
import pandas as pd
import numpy as np
from snakemake.utils import logger, validate
from snakemake.io import _load_configfile

# Determine wrapper prefix since we mix local wrappers with wrappers
# from snakemake-wrappers
SMK_WRAPPER_VERSION = "0.67.0"
SMK_WRAPPER_PREFIX_RAW = "https://github.com/snakemake/snakemake-wrappers/raw"
SMK_WRAPPER_PREFIX = f"{SMK_WRAPPER_PREFIX_RAW}/{SMK_WRAPPER_VERSION}"
WRAPPER_PREFIX = workflow.wrapper_prefix.rstrip("/")
if WRAPPER_PREFIX == SMK_WRAPPER_PREFIX_RAW:
    # Change main to version number once we start sem-versioning
    WRAPPER_PREFIX = "https://raw.githubusercontent.com/NBISweden/manticore-smk/main/workflow/wrappers"

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##################################################
# Core configuration
##################################################
include: "core/config.smk"

##### load config and sample sheets #####
configfile: "config/config.yaml"
# FIXME: Redundant at present; needed to save config
preprocess_config(config)
validate(config, schema="../schemas/config.schema.yaml")

individuals = SampleData(config["samples"].get("individual", None),
                         ignore=config["samples"].get("ignore", None))
pools = SampleData(config["samples"].get("pool", None),
                   ignore=config["samples"].get("ignore", None))
allsamples = SampleData(individuals, pools)
reads = ReadData(config["reads"].get("readfile", None),
                 ignore=config["reads"].get("ignore", None))
reads = reads.subset(samples=allsamples.unique_samples.tolist())

# FIXME: sometimes access to sample data is necessary in config
# functions (e.g. when subsetting on population) or calculating
# ploidy. Then again samples dictionary is not enormous
config["__allsamples__"] = allsamples

# Wrap config dictionary
cfg = Config(config)

##############################
## Paths - set in and get from config?
##############################
__EXTERNAL__ = Path(cfg.fs.external)
__INTERIM__ = Path(cfg.fs.interim)
__INTERIM_POOL__ = Path(cfg.fs.interim) / "pool"
__INTERIM_IND__ = Path(cfg.fs.interim) / "ind"
__METADATA__ = Path(cfg.fs.metadata)
__RAW__ = Path(cfg.fs.raw)
__RESOURCES__ = Path(cfg.fs.resources)
__REPORTS__ = Path("reports")
__RESULTS__ = Path("results")
__RESULTS_POOL__ = Path("results") / "pool"
__RESULTS_IND__ = Path("results") / "ind"

if cfg.workflow.trim:
    reads.trim(__RAW__, __INTERIM__)

##############################
## Wildcard constraints
##############################
wildcard_constraints:
    external = str(__EXTERNAL__),
    interim = str(__INTERIM__),
    interim_pool = str(__INTERIM_POOL__),
    interim_ind = str(__INTERIM_IND__),
    metadata = str(__METADATA__),
    raw = str(__RAW__),
    reports = str(__REPORTS__),
    resources = str(__RESOURCES__),
    results = str(__RESULTS__)

## File extensions
ext = {
    'bam': [".sam", ".bam", ".cram"],
    'bgz': [".bgz"],
    'fa': [".fa", ".fasta"],
    'fastq': [".fq", ".fastq"],
    'gz': [".gz", ".bgz"],
    'gatk_modes': [".g"],
    'readno': ["_1", "_2"],
    'rm': ['.rm']
}
wildcard_constraints:
    aligner = wildcards_or(["bwa"]),
    all = "all",
    analysis = wildcards_or(cfg.analysisnames),
    bam = wildcards_or(ext["bam"]),
    bamfastq = wildcards_or(ext["bam"] + ext["fastq"]),
    bgz = wildcards_or(ext["bgz"], True),
    caller = wildcards_or(cfg.workflow.variantcallers.ind + cfg.workflow.variantcallers.pool),
    callset = "rawvc",
    dot = "(|.)",
    fa = wildcards_or(ext["fa"]),
    fastq = wildcards_or(ext["fastq"]),
    filternum = "[0-9]{2}",
    filtername = wildcards_or(filter_schema.names),
    genome = os.path.splitext(os.path.basename(cfg['db']['ref']))[0],
    group = "(ind|pool)",
    gz = wildcards_or(ext["gz"], True),
    ind_vc = wildcards_or(cfg['workflow']['variantcallers']['ind']),
    kmer = "[0-9]+",
    ossep = "(|/)",
    partition = "[0-9]+",
    pool_vc = wildcards_or(cfg.workflow.variantcallers.pool),
    population = wildcards_or(allsamples.populations, empty=True),
    readno = wildcards_or(ext["readno"]),
    region = wildcards_or(cfg.workflow.regions.keys()),
    repeatmask = wildcards_or(ext["rm"], True),
    sample = wildcards_or(allsamples.unique_samples),
    sex = wildcards_or(definitions.definitions.ploidy.properties.keys()),
    statnum = "[0-9]{2}",
    statname = "(windowed_statistic|plain_statistic)",
    step_size = "[0-9]+",
    target = "[0-9]+",
    tool = wildcards_or(definitions.definitions.tool.enum)

wc = workflow._wildcard_constraints

## Store some workflow metadata
cfg["__workflow_basedir__"] = workflow.basedir
cfg["__workflow_workdir__"] = os.getcwd()
cfg["__worfklow_commit__"] = None
cfg["__worfklow_commit_link__"] = None

try:
    with cd(workflow.basedir, logger):
        commit = sp.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
        commit_short = sp.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()
        cfg["__workflow_commit__"] = commit_short
        cfg["__workflow_commit_link__"] = f"https://github.com/NBISweden/manticore-smk/commit/{commit}"
except Exception as e:
    print(e)
    raise


##################################################
## Uri parsing functions
##################################################
include: "core/uri.smk"

##################################################
# Input collection functions
##################################################
def all(wildcards):
    d = {
        'multiqc': [str(__REPORTS__ / "qc/multiqc.html")],
    }
    d.update(**all_rawvc_input(wildcards))
    d.update(**all_popoolation_input(wildcards))
    d.update(**all_popoolation2_input(wildcards))
    d.update(**all_analysisset_input(wildcards))
    #d['stats'] = all_bcftools_stats(wildcards)
    d['config'] = "results/config/manticore.config.yaml"
    return d

##############################
# qc
##############################
include: "input/qc.smk"

##############################
# Mapping
##############################
include: "input/map.smk"

##############################
# Raw variant calling
##############################
include: "input/rawvc.smk"

##############################
# WIP: BQSR
##############################

##############################
# filters
##############################
include: "input/filters.smk"

##############################
# WIP: VQSR
##############################

##############################
# Popoolation modules
##############################
include: "input/popoolation.smk"
include: "input/popoolation2.smk"

##############################
# Analysis sets
##############################
include: "input/analysisset.smk"
