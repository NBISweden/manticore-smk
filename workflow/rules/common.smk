import sys
import re
import os
import itertools
import urllib
import pandas as pd
import numpy as np
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

def _read_tsv(infile, index, schema):
    if infile is None:
        return None
    df = pd.read_csv(infile, sep="\t").set_index(index, drop=False)
    df = df.replace({np.nan: None})
    df.index.names = index
    validate(df, schema=schema)
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

individuals = _subset_by_ignore(individuals, config["samples"]["ignore"])
pools = _subset_by_ignore(pools, config["samples"]["ignore"])
reads = _subset_by_ignore(reads, config["reads"]["ignore"])
samples = individuals["SM"].tolist() + pools["SM"].tolist()
# Subset reads to samples list
reads = reads[reads.index.isin(samples, "SM")]
# Add "reads.trimmed" column to indicate reads that go into mapping
reads["reads.trimmed"] = reads["reads"]
if config["workflow"]["trim"]:
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
    'readno': ["_1", "_2"]
}

wildcard_constraints:
    bam = wildcards_or(ext["bam"]),
    bamfastq = wildcards_or(ext["bam"] + ext["fastq"]),
    bgz = wildcards_or(ext["bgz"], True),
    fa = wildcards_or(ext["fa"]),
    fastq = wildcards_or(ext["fastq"]),
    genome = config['db']['ref']['alias'],
    gz = wildcards_or(ext["gz"], True),
    kmer = "\d+",
    readno = wildcards_or(ext["readno"]),
    sample = wildcards_or(samples)

wc = workflow._wildcard_constraints

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
    return d

##############################
# qc
##############################
def trim_all(wildcards):
    if not config["workflow"]["trim"]:
        return {'reads': [], 'qc': []}
    seq = reads["reads.trimmed"]
    # List reads as paired-end or single-end
    df = reads.groupby(level=["SM", "unit"]).size().to_frame("pe")
    df = reads[reads["id"] == 1].droplevel("id").join(df["pe"])
    df["pe"] = df["pe"].map({1: "se", 2: "pe"})
    ext = rf"(_1)({wc['fastq']}{wc['gz']})$"
    df["metrics"] = df["reads"].str.replace(
        str(__RAW__), f"logs/{str(__INTERIM__)}/mapping/trim")
    df["metrics"] = df["metrics"].str.replace(
            ext, "\\2.{pe}.cutadapt_metrics.txt")
    qc = [x.format(pe=pe) for x, pe in zip(df["metrics"], df["pe"])]
    return dict(reads=seq, qc=qc)

def multiqc_all(wildcards):
    val = {
        'qc': trim_all(wildcards)['qc']  + fastqc_all(wildcards),
        'jellyfish': jellyfish_all(wildcards),
        'picard': picard_alignqc_all(wildcards)
    }
    return val

def fastqc_all(wildcards):
    # Regular input reads
    qc = [f"{str(__RESULTS__)}/qc/fastqc/{x}_fastqc.zip" for x in reads['reads']]
    if config["workflow"]["trim"]:
        qc = qc + [f"{str(__RESULTS__)}/qc/fastqc/{x}_fastqc.zip" for x in trim_all(wildcards)['reads']]
    return qc

def jellyfish_count(wildcards):
    df = reads[reads["SM"] == wildcards.sample]
    if wildcards.trimmed == ".trimmed":
        seq = df['reads'].str.replace(str(__RAW__), str(__INTERIM__ / "mapping/trim")).tolist()
    else:
        seq = df['reads'].tolist()
    return {'seq': seq}

def jellyfish_all(wildcards):
    val = []
    trim = ".trimmed" if config["workflow"]["trim"] else ""
    for kmer in config['qc']['jellyfish']['kmer']:
        hist = expand(str(__RESULTS__ / "qc/jellyfish/{sample}{trim}.{kmer}_jf.hist"),
                      sample=samples, trim=trim, kmer=kmer)
        val.extend(hist)
    return val

def picard_alignqc_all(wildcards):
    inbam = bwa_mem_all(wildcards)
    dupbam = [os.path.join(os.path.dirname(x), "dedup", os.path.basename(x)) for x in inbam]
    align_metrics = ["{}/qc/align/{}.align_metrics.txt".format(str(__RESULTS__), x) for x in inbam]
    insert_metrics = ["{}/qc/align/{}.insert_metrics.txt".format(str(__RESULTS__), x) for x in inbam]
    dup_metrics = [f"{x}.dup_metrics.txt" for x in dupbam]
    return align_metrics + insert_metrics + dup_metrics

##############################
# Mapping
##############################
bwa_ext = [".amb", ".ann", ".bwt", ".pac", ".sa"]

def bwa_mem_rg(wildcards):
    df = reads[(reads["id"] == 1) & (reads["SM"] == wildcards.sample) & (reads["unit"] == wildcards.unit)]
    rg = {
        'LB': 'lib',
        'PL': 'ILLUMINA',
        'SM': wildcards.sample,
        'PU': wildcards.unit
    }
    rglist = df.to_dict('rglist')
    assert len(rglist) == 1, "only one sample and unit should match read configuration"
    rg.update(rglist[0])
    rg['id'] = config['reads']['read_group_fmt'].format(**rglist[0])
    rgstr = r'-R "@RG\tID:{id}LB:{LB}\tPL:{PL}\tSM:{SM}\tPU:{PU}"'.format(**rg)
    return rgstr


def bwa_mem_input(wildcards):
    df = reads[(reads["SM"] == wildcards.sample) & (reads["unit"] == wildcards.unit)]
    return sorted(df['reads.trimmed'].to_list())

def bwa_mem_all(wildcards):
    df = reads[reads["id"] == 1]
    fn = str(__INTERIM__/ "map/bwa/{genome}/{{SM}}.bam")
    fn = fn.format(genome=config['db']['ref']['alias'])
    bam = [fn.format(**x) for k, x in df.iterrows()]
    return bam

def picard_merge_sam_input(wildcards):
    df = reads[(reads["SM"] == wildcards.sample) & (reads["id"] == 1)]
    fn = str(__INTERIM__/ "map/bwa/{genome}/{{SM}}/{{unit}}.bam")
    fn = fn.format(genome=config['db']['ref']['alias'])
    bam = [fn.format(**x) for k, x in df.iterrows()]
    return bam
