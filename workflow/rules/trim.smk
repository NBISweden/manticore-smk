rule all_trim:
    """Pseudo-rule to trim all input sequences"""
    input: unpack(all_trim)


rule cutadapt_pe:
    """Cutadapt: cut paired end sequences

    Cut both five- and threeprime adapter.
    """
    output:
        fastq1 = temp("{interim}/mapping/trim/{prefix}_1{fastq}{gz}"),
        fastq2 = temp("{interim}/mapping/trim/{prefix}_2{fastq}{gz}"),
        qc = "logs/{interim}/mapping/trim/{prefix}{fastq}{gz}.pe.cutadapt_metrics.txt"
    input:
        read1 = __RAW__ / "{prefix}_1{fastq}{gz}",
        read2 = __RAW__ / "{prefix}_2{fastq}{gz}"
    resources:
        runtime = lambda wildcards, attempt: attempt * config['trim']['cutadapt']['pe']['runtime']
    params:
        others = config['trim']['cutadapt']['pe']['options'],
        adapters = config['trim']['cutadapt']['pe']['adapters']
    threads: config['trim']['cutadapt']['pe']['threads']
    log: "logs/{interim}/mapping/trim/{prefix}{fastq}{gz}.log"
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/bio/cutadapt/pe"


rule cutadapt_se:
    """Cutadapt: cut single end sequences"""
    resources:
        runtime = lambda wildcards, attempt: attempt * config['trim']['cutadapt']['se']['runtime']
    params:
        " ".join([config['trim']['cutadapt']['se']['options'],
                  config['trim']['cutadapt']['se']['adapters']])
    input: __RAW__ / "{prefix}_1{fastq}{gz}",
    output:
        fastq = temp("{interim}/mapping/trim/{prefix}{fastq}{gz}"),
        qc = "logs/{interim}/mapping/trim/{prefix}{fastq}{gz}.se.cutadapt_metrics.txt"
    threads: config['trim']['cutadapt']['se']['threads']
    log: "logs/{interim}/mapping/trim/{prefix}{fastq}{gz}.log"
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/bio/cutadapt/se"
