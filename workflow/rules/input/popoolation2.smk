def all_popoolation2_input(wildcards):
    """Collact all popoolation2 targets"""
    val = {}
    if len(pools) == 0:
        return val
    val.update(**all_popoolation2_raw_input(wildcards))
    return val


def all_popoolation2_raw_input(wildcards):
    """Generate all raw popoolation2 targets"""
    if len(pools) == 0:
        return {}
    pfx = str(__RESULTS_POOL__ / "raw/popoolation2/{sex}.{region}.{partition}.sync.gz")
    val = {'popoolation2.raw': []}
    for region in config["workflow"]["regions"].keys():
        npart = config["workflow"]["regions"][region]["npart"]
        for sex in config["workflow"]["regions"][region]["ploidy"].keys():
            val['popoolation2.raw'].extend(expand(pfx, region=region, sex=sex,
                                                  partition=list(range(npart))))
    return val


def popoolation2_samtools_mpileup_input(wildcards):
    fn = str(__INTERIM__ / "map/bwa/{sample}.bam")
    df = pools
    if wildcards.sex != "common":
        df = pools[pools.sex.isin([wildcards.sex])]
    bam = expand(fn, sample = df.SM.tolist())
    bai = [f"{x}.bai" for x in bam]
    ref = config['db']['ref']
    targets = os.path.join(
        os.path.dirname(ref), "popoolation", f"{wildcards.region}.{wildcards.target}.bed")
    return {'bam': bam, 'bai': bai, 'targets': targets}


def popoolation2_gather_parallel_results_input(wildcards):
    npart = config['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    analysis = expand("{{interim_pool}}/popoolation2/{{region}}/{{prefix}}{{repeatmask}}.{{filters}}{partition}.{{analysis}}.{{suffix}}",
                      partition = list(range(npart)), **wildcards)
    return {'analysis': analysis}


def popoolation2_filter_pileup_by_gtf_input(wildcards):
    """Generate indels input file from popoolation2, all samples"""
    d = get_filter_input(wildcards)
    if "gtf" in d.keys():
        return d["gtf"]
    return "{results}/pool/raw/popoolation2.indels/common.{region}.{target}.mpileup.gz.indels.gtf".format(**wildcards)
