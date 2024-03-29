def all_map_input(wildcards):
    targets = (
        all_fastqc(wildcards)
        + all_jellyfish(wildcards)
        + list(
            set(
                all_map_sample_targets(wildcards)
                + all_map_sample_dedup_targets(wildcards)
            )
        )
    )
    return targets


bwa_ext = [".amb", ".ann", ".bwt", ".pac", ".sa"]


def bwa_mem_rg(wildcards):
    df = reads.subset(id=1, samples=wildcards.sample, unit=wildcards.unit).data
    rg = {"LB": "lib", "PL": "ILLUMINA", "SM": wildcards.sample, "PU": wildcards.unit}
    rglist = df.to_dict("rglist")
    assert len(rglist) == 1, "only one sample and unit should match read configuration"
    rg.update(rglist[0])
    rg["id"] = cfg["reads"]["read_group_fmt"].format(**rglist[0])
    rgstr = r'-R "@RG\tID:{id}LB:{LB}\tPL:{PL}\tSM:{SM}\tPU:{PU}"'.format(**rg)
    return rgstr


def bwa_mem_index(wildcards):
    """Retrieve bma mem index"""
    return f"{wildcards.interim}/map/bwa/index/{cfg.genome}"


def bwa_mem_index_ext(wildcards):
    """Retrieve bwa mem index with extensions"""
    index = bwa_mem_index(wildcards)
    return expand("{index}{ext}", index=index, ext=bwa_ext)


def map_sample_unit_input(wildcards):
    """Retrieve input to a mapping job for a given sample and unit"""
    df = reads.subset(SM=wildcards.sample, unit=wildcards.unit)
    return sorted(df.reads.to_list())


def map_sample_target(wildcards):
    """Retrieve mapping sample target name for a given alignment program"""
    df = reads.subset(id=1, samples=wildcards.sample)
    fn = str(__INTERIM__ / f"map/bwa/{{SM}}.bam")
    bam = [fn.format(SM=x) for x in df.unique_samples]
    return bam


def map_sample_target_bai(wildcards):
    """Retrieve mapping sample index target name for a given alignment program"""
    bam = map_sample_target(wildcards)
    return [f"{x}.bai" for x in bam]


def map_dedup_sample_target(wildcards):
    """Return *unique* bam dedup sample targets"""
    df = reads.subset(id=1, samples=wildcards.sample)
    fn = str(__INTERIM__ / "map/bwa/dedup/{SM}.bam")
    bam = [fn.format(SM=sm) for sm in df.unique_samples.tolist()]
    return bam


def all_map_sample_targets(wildcards):
    """All merged map sample targets for a given alignment program"""
    df = reads.subset(id=1)
    fn = str(__INTERIM__ / "map/bwa/{SM}.bam")
    bam = [fn.format(SM=sm) for sm in df.unique_samples.tolist()]
    return bam


def all_map_sample_dedup_targets(wildcards):
    """All merged and deduplicated map sample targets for a given alignment program.

    Does not apply to pools.
    """
    df = reads.subset(id=1, samples=individuals.unique_samples.tolist())
    fn = str(__INTERIM__ / "map/bwa/dedup/{SM}.bam")
    bam = [fn.format(SM=sm) for sm in df.unique_samples.tolist()]
    return bam


def map_picard_merge_sam_input(wildcards):
    """Generate picard merge input file targets

    The map_picard_merge_sam rule merges mapped units to a
    sample-level bam file.

    """
    df = reads.subset(samples=wildcards.sample, id=1)
    fn = str(__INTERIM__ / "map/bwa/{SM}/{unit}.bam")
    bam = list(set([fn.format(**x) for k, x in df.data.iterrows()]))
    return bam
