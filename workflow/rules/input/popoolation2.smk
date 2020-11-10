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
