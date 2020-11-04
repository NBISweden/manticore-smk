def gatkhc_targets_input(wildcards):
    ref = config['db']['ref']
    fai = f"{ref}.fai"
    targets = os.path.join(
        os.path.dirname(ref), "gatkhc", f"{wildcards.region}.{wildcards.target}.bed")
    bam = bwa_mem_sample(wildcards)
    bai = [f"{x}.bai" for x in bam]
    faext = wildcards_or(ext["fa"])
    d = re.sub(faext, ".dict", ref)
    return {'ref': ref, 'fai': fai, 'targets': targets,
            'bam': bam, 'bai': bai, 'dict': d}
