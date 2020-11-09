import sys
import re
import os
import itertools
import urllib
import pandas as pd
import numpy as np
from snakemake.io import _load_configfile
from snakemake.utils import logger, validate

# Determine wrapper prefix since we mix local wrappers with wrappers
# from snakemake-wrappers
SMK_WRAPPER_VERSION = "0.67.0"
SMK_WRAPPER_PREFIX_RAW = "https://github.com/snakemake/snakemake-wrappers/raw"
SMK_WRAPPER_PREFIX = f"{SMK_WRAPPER_PREFIX_RAW}/{SMK_WRAPPER_VERSION}"
WRAPPER_PREFIX = workflow.wrapper_prefix
if WRAPPER_PREFIX == SMK_WRAPPER_PREFIX_RAW:
    # Change main to version number once we start sem-versioning
    WRAPPER_PREFIX = "https://raw.githubusercontent.com/nbis/manticore-smk/main/workflow/wrappers"


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

##############################
## Paths
##############################
__EXTERNAL__ = Path(config["fs"]["external"])
__INTERIM__ = Path(config["fs"]["interim"])
__METADATA__ = Path(config["fs"]["metadata"])
__RAW__ = Path(config["fs"]["raw"])
__REPORTS__ = Path("reports")
__RESOURCES__ = Path(config["fs"]["resources"])
__RESULTS__ = Path("results")

def _read_tsv(infile, index, schemafile):
    if infile is None:
        schema = _load_configfile(os.path.join(workflow.current_basedir, schemafile))
        cols = schema["required"]
        df = pd.DataFrame({k:[] for k in cols})
        df = df.set_index(index, drop=False)
    else:
        df = pd.read_csv(infile, sep="\t").set_index(index, drop=False)
        df = df.replace({np.nan: None})
    df.index.names = index
    validate(df, schema=schemafile)
    return df

individuals = _read_tsv(config["samples"]["individual"], ["SM"],
                        "../schemas/samples.ind.schema.yaml")
pools = _read_tsv(config["samples"]["pool"], ["SM"],
                  "../schemas/samples.pool.schema.yaml")
datasources = _read_tsv(config["datasources"], ["data"],
                        "../schemas/datasources.schema.yaml")
reads = _read_tsv(config["reads"]["readfile"], ["SM", "unit", "id"],
                  "../schemas/reads.schema.yaml")


# Subset by ignore
def _subset_by_ignore(df, ignore):
    if ignore is None:
        return df
    i = df.index.isin(ignore)
    if any(i.tolist()):
        logger.info(f"Removing entries {df[i]}")
        df = df[~i]
    return df

individuals = _subset_by_ignore(individuals, config["samples"].get("ignore", []))
pools = _subset_by_ignore(pools, config["samples"].get("ignore", []))
reads = _subset_by_ignore(reads, config["reads"].get("ignore", []))
samples = individuals["SM"].tolist() + pools["SM"].tolist()
# Subset reads to samples list
reads = reads[reads.index.isin(samples, "SM")]
# Add "reads.trimmed" column to indicate reads that go into mapping
reads["reads.trimmed"] = reads["reads"]
if config["workflow"]["trim"] and len(reads) > 0:
    reads["reads.trimmed"] = reads["reads"].str.replace(str(__RAW__), str(__INTERIM__ / "mapping/trim"))

# Save current base dir for later validation in functions
BASEDIR = workflow.current_basedir


##############################
## Wildcard constraints
##############################
wildcard_constraints:
    external = str(__EXTERNAL__),
    interim = str(__INTERIM__),
    metadata = str(__METADATA__),
    raw = str(__RAW__),
    reports = str(__REPORTS__),
    resources = str(__RESOURCES__),
    results = str(__RESULTS__)

def wildcards_or(items, empty=False):
    if empty:
        items = [""] + items
    return f'({"|".join(items)})'

## File extensions
ext = {
    'bam': [".sam", ".bam", ".cram"],
    'bgz': [".bgz"],
    'fa': [".fa", ".fasta"],
    'fastq': [".fq", ".fastq"],
    'gz': [".gz", ".bgz"],
    'gatk_modes': [".g"],
    'readno': ["_1", "_2"]
}
wildcard_constraints:
    bam = wildcards_or(ext["bam"]),
    bamfastq = wildcards_or(ext["bam"] + ext["fastq"]),
    bgz = wildcards_or(ext["bgz"], True),
    fa = wildcards_or(ext["fa"]),
    fastq = wildcards_or(ext["fastq"]),
    genome = os.path.splitext(os.path.basename(config['db']['ref']))[0],
    gz = wildcards_or(ext["gz"], True),
    kmer = "[0-9]+",
    partitions = "[0-9]+",
    readno = wildcards_or(ext["readno"]),
    region = wildcards_or(config['workflow']['regions'].keys()),
    sample = wildcards_or(samples),
    ind_vc = wildcards_or(config['workflow']['variantcallers']['ind'])

wc = workflow._wildcard_constraints

##################################################
# Core configuration
##################################################
include: "core/config.smk"

##################################################
## Uri parsing functions
##################################################
def get_uri_scheme(uri):
    return urllib.parse.urlparse(uri).scheme

def get_uri_netloc(uri):
    return urllib.parse.urlparse(uri).netloc

def parse_uri(uri):
    """Parse uri and return snakemake target"""
    scheme = get_uri_scheme(uri)
    uri = re.sub(f"{scheme}://", "", uri)
    if not scheme in ['', 'rsync', 'file', 'sftp']:
        logger.error(f"scheme '{scheme}' not allowed: use one of '', 'file', 'rsync' or 'sftp'")
        sys.exit(1)
    if scheme in ['', 'file'] and not uri.startswith("/"):
        uri = os.path.normpath(os.path.abspath(uri))
    if scheme == 'sftp':
        try:
            from snakemake.remote.SFTP import RemoteProvider
            SFTP = RemoteProvider()
            uri = SFTP.remote(uri)
        except WorkflowError as e:
            logger.error(e)
    return uri

def manticore_get_external_input(uri):
    if get_uri_scheme(uri) == "sftp":
        return parse_uri(uri)
    netloc = get_uri_netloc(uri)
    uri = parse_uri(uri)
    if re.search("[:@]", netloc):
        return []
    return uri

##################################################
# Input collection functions
##################################################
def all(wildcards):
    d = {
        'multiqc': [str(__REPORTS__ / "qc/multiqc.html")],
    }
    d.update(**all_rawvc(wildcards))
    return d

##############################
# qc
##############################
include: "input/qc.smk"

##############################
# Mapping
##############################
include: "input/mapping.smk"

##############################
# Raw variant calling
##############################
include: "input/rawvc.smk"
