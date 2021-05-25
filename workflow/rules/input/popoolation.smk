def all_popoolation_input(wildcards):
    """Collect all popoolation targets"""
    val = {}
    if len(pools.data) == 0:
        return val
    val.update(**all_popoolation_raw_input(wildcards))
    return val


def all_popoolation_raw_input(wildcards):
    """Generate popoolation raw pileup files"""
    if len(pools.data) == 0:
        return {}
    pfx = str(__RESULTS_POOL__ / "raw/popoolation/{SM}.{region}.{partition}.pileup.gz")
    val = {'popoolation.raw': []}
    for region in cfg["workflow"]["regions"].keys():
        npart = cfg["workflow"]["regions"][region]["npart"]
        val['popoolation.raw'].extend(expand(pfx, SM=pools.samples, region=region, partition=list(range(npart))))
    return val


def popoolation_gather_parallel_results_input(wildcards):
    """Generate input files for scattered popoolation jobs"""
    npart = cfg['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    analysis = expand("{{results}}/{{group}}/analysis/{{analysis}}/s{{itemnum}}_{{statname}}_{{tool}}/{{sample}}.{{region}}.w{{window_size}}.s{{step_size}}.{{measure}}.{partition}.txt.gz",
                      partition = list(range(npart)), **wildcards)
    return {'analysis': analysis}


def popoolation_filter_pileup_by_gtf_input(wildcards):
    """Generate gtf input file from popoolation2, all samples"""
    return "{results}/pool/raw/popoolation2.indels/common.{region}.{target}.mpileup.gz.indels.gtf".format(**wildcards)


def popoolation_samtools_filter_mpileup_input(wildcards):
    ref = cfg['db']['ref']
    targets = os.path.join(
        os.path.dirname(ref), "popoolation", f"{wildcards.region}.{wildcards.target}.bed")
    bam = str(__INTERIM__ / f"map/bwa/{wildcards.sample}.bam")
    return {'bam': bam, 'targets': targets}
