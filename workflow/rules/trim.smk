rule all_trim:
    """Pseudo-rule to trim all input sequences"""
    input: lambda wildcards: all_fastqc(wildcards) + all_jellyfish(wildcards) + all_trim(wildcards)['qc']


rule trim_cutadapt_pe:
    """Cutadapt: cut paired end sequences

    Cut both five- and threeprime adapter.
    """
    output:
        fastq1 = temp("{interim}/map/trim/{prefix}_1{fastq}{gz}"),
        fastq2 = temp("{interim}/map/trim/{prefix}_2{fastq}{gz}"),
        qc = "logs/{interim}/map/trim/{prefix}{fastq}{gz}.pe.cutadapt_metrics.txt"
    input:
        read1 = __RAW__ / "{prefix}_1{fastq}{gz}",
        read2 = __RAW__ / "{prefix}_2{fastq}{gz}"
    resources:
        runtime = lambda wildcards, attempt: cfg.ruleconf("trim_cutadapt_pe", attempt).resources("runtime")
    params:
        others = cfg.ruleconf("trim_cutadapt_pe").params("options"),
        adapters = cfg.ruleconf("trim_cutadapt_pe").params("adapters")
    threads: cfg.ruleconf("trim_cutadapt_pe").threads
    log: "logs/{interim}/map/trim/{prefix}{fastq}{gz}.log"
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/bio/cutadapt/pe"


rule trim_cutadapt_se:
    """Cutadapt: cut single end sequences"""
    output:
        fastq = temp("{interim}/map/trim/{prefix}{fastq}{gz}"),
        qc = "logs/{interim}/map/trim/{prefix}{fastq}{gz}.se.cutadapt_metrics.txt"
    input: __RAW__ / "{prefix}_1{fastq}{gz}"
    resources:
        runtime = lambda wildcards, attempt: cfg.ruleconf("trim_cutadapt_se", attempt).resources("runtime")
    params:
        " ".join([cfg.ruleconf("trim_cutadapt_se").params("options"),
                  cfg.ruleconf("trim_cutadapt_se").params("adapters")])
    threads: cfg.ruleconf("trim_cutadapt_se").threads
    log: "logs/{interim}/map/trim/{prefix}{fastq}{gz}.log"
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/bio/cutadapt/se"
