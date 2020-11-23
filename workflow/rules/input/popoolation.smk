def all_popoolation_input(wildcards):
    """Collect all popoolation targets"""
    val = {}
    if len(pools) == 0:
        return val
    val.update(**all_popoolation_raw_input(wildcards))
    val.update(**all_popoolation_analysis_input(wildcards))
    return val


def all_popoolation_raw_input(wildcards):
    """Generate popoolation raw pileup files"""
    if len(pools) == 0:
        return {}
    pfx = str(__RESULTS_POOL__ / "raw/popoolation/{SM}.{region}.{partition}.pileup.gz")
    val = {'popoolation.raw': []}
    for region in config["workflow"]["regions"].keys():
        npart = config["workflow"]["regions"][region]["npart"]
        val['popoolation.raw'].extend(expand(pfx, SM=pools.SM, region=region, partition=list(range(npart))))
    return val


def all_popoolation_analysis_input(wildcards):
    """Generate popoolation analysis stats and plot results"""
    if len(pools) == 0:
        return {}
    val = {'popoolation.analysis': []}
    for k in config.keys():
        if not k.startswith("analysis/"):
            continue
        pfx = __RESULTS_POOL__ / f"{k}/{{stat}}/{{sample}}.{{region}}.w{{window_size}}.s{{step_size}}.{{measure}}.txt.gz"
        regions = list(config["workflow"]["regions"].keys())
        statistics = config[k].get("statistics", [])
        if len(statistics) == 0:
            continue
        for k, stat in statistics.items():
            if stat['engine'] != 'popoolation':
                continue
            wsize = [10000]
            res = expand(pfx, sample=pools["SM"], stat=k, region=regions, window_size=wsize, step_size=wsize, measure=stat['statistic'])
            val['popoolation.analysis'].extend(res)

    return val





def popoolation_gather_parallel_results_input(wildcards):
    """Generate input files for scattered popoolation jobs"""
    npart = config['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    analysis = expand("{{results}}/pool/analysis/{{analysis}}/{{statname}}/{{sample}}.{{region}}.w{{window_size}}.s{{step_size}}.{{measure}}.{partition}.txt.gz",
                      partition = list(range(npart)), **wildcards)
    return {'analysis': analysis}


def popoolation_filter_pileup_by_gtf_input(wildcards):
    return "{results}/pool/raw/popoolation.indels/common.{region}.{target}.mpileup.gz.indels.gtf".format(**wildcards)


def popoolation_samtools_filter_mpileup_input(wildcards):
    ref = config['db']['ref']
    targets = os.path.join(
        os.path.dirname(ref), "popoolation", f"{wildcards.region}.{wildcards.target}.bed")
    bam = str(__INTERIM__ / f"map/bwa/{wildcards.sample}.bam")
    return {'bam': bam, 'targets': targets}
