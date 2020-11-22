import sys
import re
import os
import itertools
import urllib
import contextlib
import subprocess as sp
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
__INTERIM_POOL__ = Path(config["fs"]["interim"]) / "pool"
__INTERIM_IND__ = Path(config["fs"]["interim"]) / "ind"
__METADATA__ = Path(config["fs"]["metadata"])
__RAW__ = Path(config["fs"]["raw"])
__REPORTS__ = Path("reports")
__RESOURCES__ = Path(config["fs"]["resources"])
__RESULTS__ = Path("results")
__RESULTS_POOL__ = Path("results") / "pool"
__RESULTS_IND__ = Path("results") / "ind"

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
    ## Check for duplicate rows
    if any(df.duplicated()):
        logger.error(f"duplicated rows in {infile}:")
        logger.error(f"  {df[df.duplicated()]}")
        sys.exit(1)
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
        logger.info(f"Ignoring entries in this analysis: {df[i]}")
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
    reads["reads.trimmed"] = reads["reads"].str.replace(str(__RAW__), str(__INTERIM__ / "map/trim"))

# Save current base dir for later validation in functions
BASEDIR = workflow.current_basedir


##################################################
# Core configuration
##################################################
include: "core/config.smk"

##############################
## Wildcard constraints
##############################

##### load schema definitions #####
def load_schema(schema):
    import inspect
    import yaml
    if not os.path.isabs(schema):
        frame = inspect.currentframe().f_back
    if "workflow" in frame.f_globals:
        workflow = frame.f_globals["workflow"]
        schema = os.path.join(workflow.current_basedir, schema)
    with open(schema) as fh:
        data = yaml.safe_load(fh)
    return data

definitions = load_schema("../schemas/definitions.schema.yaml")

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
    'readno': ["_1", "_2"],
    'rm': ['.rm']
}
wildcard_constraints:
    aligner = wildcards_or(["bwa"]),
    all = "all",
    analysis = wildcards_or(get_analysisnames()),
    bam = wildcards_or(ext["bam"]),
    bamfastq = wildcards_or(ext["bam"] + ext["fastq"]),
    bgz = wildcards_or(ext["bgz"], True),
    caller = wildcards_or(config['workflow']['variantcallers']['ind'] + config['workflow']['variantcallers']['pool']),
    callset = "rawvc",
    fa = wildcards_or(ext["fa"]),
    fastq = wildcards_or(ext["fastq"]),
    filtername = wildcards_or(get_filternames()),
    genome = os.path.splitext(os.path.basename(config['db']['ref']))[0],
    group = "(ind|pool)",
    gz = wildcards_or(ext["gz"], True),
    ind_vc = wildcards_or(config['workflow']['variantcallers']['ind']),
    kmer = "[0-9]+",
    ossep = "(|/)",
    partition = "[0-9]+",
    pool_vc = wildcards_or(config['workflow']['variantcallers']['pool']),
    readno = wildcards_or(ext["readno"]),
    region = wildcards_or(config['workflow']['regions'].keys()),
    repeatmask = wildcards_or(ext["rm"], True),
    sample = wildcards_or(samples),
    sex = wildcards_or(list(definitions["definitions"]["ploidy"]["properties"].keys())),
    step_size = "[0-9]+",
    target = "[0-9]+",

wc = workflow._wildcard_constraints

## Store some workflow metadata
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_workdir__"] = os.getcwd()
config["__worfklow_commit__"] = None

try:
    with cd(workflow.basedir, logger):
        logger.info(f"cd to {workflow.basedir}")
        out = sp.check_output(["git", "rev-parse", "--short", "HEAD"])
        config["__worfklow_commit__"] = out.decode().strip()
except Exception as e:
    print(e)
    raise

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
    d.update(**all_rawvc_input(wildcards))
    d.update(**all_popoolation_input(wildcards))
    d.update(**all_popoolation2_input(wildcards))
    #d['stats'] = all_bcftools_stats(wildcards)
    d['config'] = "config/manticore.config.yaml"
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
# Variant filtering with hard filters
##############################
include: "input/variantfilter.smk"

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
