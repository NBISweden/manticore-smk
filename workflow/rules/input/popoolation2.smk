def popoolation2_samtools_mpileup_input(wildcards):
    fn = str(__INTERIM__ / "map/bwa/{sample}.bam")
    df = pools
    if "sex" in dict(wildcards):
        if wildcards.sex != "":
            df = pools[pools.sex.isin([wildcards.sex])]
    bam = expand(fn, sample = df.SM.tolist())
    bai = [f"{x}.bai" for x in bam]
    ref = config['db']['ref']
    targets = os.path.join(
        os.path.dirname(ref), "popoolation", f"{wildcards.region}{wildcards.repeatmask}.{wildcards.target}.bed")
    return {'bam': bam, 'bai': bai, 'targets': targets}


def popoolation2_gather_parallel_results_input(wildcards):
    npart = config['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    analysis = expand("{{interim_pool}}/popoolation2/{{region}}/{{prefix}}{{repeatmask}}.{{filters}}{partition}.{{analysis}}.{{suffix}}",
                      partition = list(range(npart)), **wildcards)
    return {'analysis': analysis}
