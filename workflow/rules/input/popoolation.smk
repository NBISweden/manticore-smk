def popoolation_filter_pileup_by_gtf_input(wildcards):
    return "{interim}/pool/popoolation2/{region}/all{repeatmask}.mpileup.gz.indels.gtf".format(**wildcards)


def popoolation_samtools_filter_mpileup_input(wildcards):
    ref = config['db']['ref']
    targets = os.path.join(
        os.path.dirname(ref), "popoolation", f"{wildcards.region}{wildcards.repeatmask}.{wildcards.target}.bed")
    bam = str(__INTERIM__ / f"map/bwa/{wildcards.sample}.bam")
    return {'bam': bam, 'targets': targets}
