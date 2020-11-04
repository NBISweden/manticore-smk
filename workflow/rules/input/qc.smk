def trim_all(wildcards):
    if not config["workflow"]["trim"]:
        return {'reads': [], 'qc': []}
    seq = reads["reads.trimmed"]
    # List reads as paired-end or single-end
    df = reads.groupby(level=["SM", "unit"]).size().to_frame("pe")
    df = reads[reads["id"] == 1].droplevel("id").join(df["pe"])
    df["pe"] = df["pe"].map({1: "se", 2: "pe"})
    ext = rf"(_1)({wc['fastq']}{wc['gz']})$"
    df["metrics"] = df["reads"].str.replace(
        str(__RAW__), f"logs/{str(__INTERIM__)}/mapping/trim")
    df["metrics"] = df["metrics"].str.replace(
            ext, "\\2.{pe}.cutadapt_metrics.txt")
    qc = [x.format(pe=pe) for x, pe in zip(df["metrics"], df["pe"])]
    return dict(reads=seq, qc=qc)

def multiqc_all(wildcards):
    val = {
        'qc': trim_all(wildcards)['qc'] + fastqc_all(wildcards),
        'jellyfish': jellyfish_all(wildcards),
        'picard': picard_alignqc_all(wildcards)
    }
    return val

def fastqc_all(wildcards):
    # Regular input reads
    qc = [f"{str(__RESULTS__)}/qc/fastqc/{x}_fastqc.zip" for x in reads['reads']]
    if config["workflow"]["trim"]:
        qc = qc + [f"{str(__RESULTS__)}/qc/fastqc/{x}_fastqc.zip" for x in trim_all(wildcards)['reads']]
    return qc

def jellyfish_count(wildcards):
    df = reads[reads["SM"] == wildcards.sample]
    if wildcards.trimmed == ".trimmed":
        seq = df['reads'].str.replace(str(__RAW__), str(__INTERIM__ / "mapping/trim")).tolist()
    else:
        seq = df['reads'].tolist()
    return {'seq': seq}

def jellyfish_all(wildcards):
    val = []
    trim = ".trimmed" if config["workflow"]["trim"] else ""
    for kmer in config['qc']['jellyfish']['kmer']:
        hist = expand(str(__RESULTS__ / "qc/jellyfish/{sample}{trim}.{kmer}_jf.hist"),
                      sample=samples, trim=trim, kmer=kmer)
        val.extend(hist)
    return val

def picard_alignqc_all(wildcards):
    inbam = bwa_mem_all(wildcards)
    dupbam = [os.path.join(os.path.dirname(x), "dedup", os.path.basename(x)) for x in inbam]
    align_metrics = ["{}/qc/align/{}.align_metrics.txt".format(str(__RESULTS__), x) for x in inbam]
    insert_metrics = ["{}/qc/align/{}.insert_metrics.txt".format(str(__RESULTS__), x) for x in inbam]
    dup_metrics = [f"{x}.dup_metrics.txt" for x in dupbam]
    return align_metrics + insert_metrics + dup_metrics
