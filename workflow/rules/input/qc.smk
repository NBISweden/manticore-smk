def all_multiqc(wildcards):
    val = {
        'qc': all_trim(wildcards)['qc'] + all_fastqc(wildcards),
        'jellyfish': all_jellyfish(wildcards),
        'picard': all_picard_alignqc(wildcards),
        'qualimap': all_qualimap_bamqc(wildcards),
        'bcftools': all_bcftools_stats(wildcards),
        'sambamba': all_sambamba_depth(wildcards)  # Should go to report.smk as we need to make separate plot of this
    }
    return val


def all_trim(wildcards):
    if not cfg.workflow["trim"]:
        return {'reads': [], 'qc': []}
    seq = reads.data["reads.trimmed"].tolist()
    df = reads.read_pairs
    ext = rf"(_1)({wc['fastq']}{wc['gz']})$"
    df["metrics"] = df["reads"].str.replace(
        str(__RAW__), f"logs/{str(__INTERIM__)}/map/trim")
    df["metrics"] = df["metrics"].str.replace(
            ext, "\\2.{pe}.cutadapt_metrics.txt")
    qc = [x.format(pe=pe) for x, pe in zip(df["metrics"], df["pe"])]
    return dict(reads=seq, qc=qc)


def all_fastqc(wildcards):
    if 'fastqc' not in cfg.workflow['qc']:
        return []
    # Regular input reads
    qc = [f"{str(__RESULTS__)}/qc/fastqc/{x}_fastqc.zip" for x in reads.data['reads'].tolist()]
    if cfg.workflow["trim"]:
        qc = qc + [f"{str(__RESULTS__)}/qc/fastqc/{x}_fastqc.zip" for x in all_trim(wildcards)['reads']]
    return qc

def jellyfish_count(wildcards):
    df = reads.subset(samples=wildcards.sample)
    if wildcards.trimmed == ".trimmed":
        seq = df.data['reads'].str.replace(str(__RAW__), str(__INTERIM__ / "map/trim")).tolist()
    else:
        seq = df.data['reads'].tolist()
    return {'seq': seq}

def all_jellyfish(wildcards):
    if 'jellyfish' not in cfg.workflow['qc']:
        return []
    val = []
    trim = ".trimmed" if cfg.workflow["trim"] else ""
    for kmer in cfg.ruleconf("qc_jellyfish").params("kmer"):
        hist = expand(str(__RESULTS__ / "qc/jellyfish/{sample}{trim}.{kmer}_jf.hist"),
                      sample=allsamples.samples, trim=trim, kmer=kmer)
        val.extend(hist)
    return val

def all_picard_alignqc(wildcards):
    if 'picard' not in cfg.workflow['qc']:
        return []
    inbam = all_map_sample_targets(wildcards)
    dupbam = all_map_sample_dedup_targets(wildcards)
    align_metrics = ["{}/qc/align/{}.align_metrics.txt".format(str(__RESULTS__), x) for x in inbam]
    insert_metrics = ["{}/qc/align/{}.insert_metrics.txt".format(str(__RESULTS__), x) for x in inbam]
    dup_metrics = [f"{x}.dup_metrics.txt" for x in dupbam]
    return align_metrics + insert_metrics + dup_metrics


def all_qualimap_bamqc(wildcards):
    if 'qualimap' not in cfg.workflow['qc']:
        return []
    inbam = all_map_sample_targets(wildcards)
    bamdf = {os.path.splitext(os.path.basename(x))[0]:os.path.basename(x) for x in inbam}
    df = {x["SM"]:x["pe"] for k, x in reads.read_pairs.iterrows()}
    fn = str(__RESULTS__/ "qc/qualimap/{bam}.{pe}.qualimap/genome_results.txt")
    results = [fn.format(bam=bamdf[SM], pe=df[SM]) for SM in bamdf.keys()]
    return results


def all_sambamba_depth(wildcards):
    """Run sambamba depth on 100bp windows"""
    if 'sambamba' not in cfg.workflow['qc']:
        return []
    inbam = all_map_sample_targets(wildcards)
    bampfx = [re.sub(wc["bam"], "", x) for x in inbam]
    fn = str(__RESULTS__/ f"qc/sambamba/{{SM}}.depth.w100.bed.gz")
    results = [fn.format(SM=x) for x in allsamples.unique_samples]
    return results


def all_bcftools_stats(wildcards):
    if 'bcftools' not in cfg['workflow']['qc']:
        return []
    invcf = all_rawvc_input(wildcards)['rawvc']
    bn = [re.sub(rf"^{str(__RESULTS__)}", str(__RESULTS__ / "qc/variants"), x) for x in invcf]
    stats = [f"{x}.stats" for x in bn]
    return stats
