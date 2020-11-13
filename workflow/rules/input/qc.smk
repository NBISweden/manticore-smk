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
    if not config["workflow"]["trim"]:
        return {'reads': [], 'qc': []}
    seq = reads["reads.trimmed"]
    df = read_pairs_dataframe()
    ext = rf"(_1)({wc['fastq']}{wc['gz']})$"
    df["metrics"] = df["reads"].str.replace(
        str(__RAW__), f"logs/{str(__INTERIM__)}/map/trim")
    df["metrics"] = df["metrics"].str.replace(
            ext, "\\2.{pe}.cutadapt_metrics.txt")
    qc = [x.format(pe=pe) for x, pe in zip(df["metrics"], df["pe"])]
    return dict(reads=seq, qc=qc)


def all_fastqc(wildcards):
    if 'fastqc' not in config['workflow']['qc']:
        return []
    # Regular input reads
    qc = [f"{str(__RESULTS__)}/qc/fastqc/{x}_fastqc.zip" for x in reads['reads']]
    if config["workflow"]["trim"]:
        qc = qc + [f"{str(__RESULTS__)}/qc/fastqc/{x}_fastqc.zip" for x in all_trim(wildcards)['reads']]
    return qc

def jellyfish_count(wildcards):
    df = reads[reads["SM"] == wildcards.sample]
    if wildcards.trimmed == ".trimmed":
        seq = df['reads'].str.replace(str(__RAW__), str(__INTERIM__ / "map/trim")).tolist()
    else:
        seq = df['reads'].tolist()
    return {'seq': seq}

def all_jellyfish(wildcards):
    if 'jellyfish' not in config['workflow']['qc']:
        return []
    val = []
    trim = ".trimmed" if config["workflow"]["trim"] else ""
    for kmer in config['resources']['qc_jellyfish']['kmer']:
        hist = expand(str(__RESULTS__ / "qc/jellyfish/{sample}{trim}.{kmer}_jf.hist"),
                      sample=samples, trim=trim, kmer=kmer)
        val.extend(hist)
    return val

def all_picard_alignqc(wildcards):
    if 'picard' not in config['workflow']['qc']:
        return []
    inbam = all_map_sample_targets(wildcards)
    dupbam = all_map_sample_dedup_targets(wildcards)
    align_metrics = ["{}/qc/align/{}.align_metrics.txt".format(str(__RESULTS__), x) for x in inbam]
    insert_metrics = ["{}/qc/align/{}.insert_metrics.txt".format(str(__RESULTS__), x) for x in inbam]
    dup_metrics = [f"{x}.dup_metrics.txt" for x in dupbam]
    return align_metrics + insert_metrics + dup_metrics


def all_qualimap_bamqc(wildcards):
    if 'qualimap' not in config['workflow']['qc']:
        return []
    inbam = all_map_sample_targets(wildcards)
    bamdf = {os.path.splitext(os.path.basename(x))[0]:os.path.basename(x) for x in inbam}
    df = {x["SM"]:x["pe"] for k, x in read_pairs_dataframe().iterrows()}
    fn = str(__RESULTS__/ "qc/qualimap/{bam}.{pe}.qualimap/genome_results.txt")
    results = [fn.format(bam=bamdf[SM], pe=df[SM]) for SM in bamdf.keys()]
    return results


def all_sambamba_depth(wildcards):
    if 'sambamba' not in config['workflow']['qc']:
        return []
    inbam = all_map_sample_targets(wildcards)
    bampfx = [re.sub(wc["bam"], "", x) for x in inbam]
    fn = str(__RESULTS__/ f"qc/sambamba/{{SM}}.depth.w100.bed.gz")
    results = [fn.format(SM=x) for x in samples]
    return results


def all_bcftools_stats(wildcards):
    if 'bcftools' not in config['workflow']['qc']:
        return []
    invcf = all_variantfilter(wildcards)['variantfilter'] + all_rawvc(wildcards)['rawvc']
    bn = [re.sub(rf"^{str(__RESULTS__)}", str(__RESULTS__ / "qc/variantstats"), x) for x in invcf]
    stats = [f"{x}.stats" for x in bn]
    return stats
