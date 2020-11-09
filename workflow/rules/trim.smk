rule all_trim:
    """Pseudo-rule to trim all input sequences"""
    input: unpack(all_trim)


rule trim_cutadapt_pe:
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
        runtime = lambda wildcards, attempt: resources("trim_cutadapt_pe", "runtime", attempt)
    params:
        others = get_params("trim_cutadapt_pe", "options"),
        adapters = get_params("trim_cutadapt_pe", "adapters")
    threads: get_params("trim_cutadapt_pe", 'threads')
    log: "logs/{interim}/mapping/trim/{prefix}{fastq}{gz}.log"
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/bio/cutadapt/pe"


rule trim_cutadapt_se:
    """Cutadapt: cut single end sequences"""
    output:
        fastq = temp("{interim}/mapping/trim/{prefix}{fastq}{gz}"),
        qc = "logs/{interim}/mapping/trim/{prefix}{fastq}{gz}.se.cutadapt_metrics.txt"
    input: __RAW__ / "{prefix}_1{fastq}{gz}"
    resources:
        runtime = lambda wildcards, attempt: resources("trim_cutadapt_se", "runtime", attempt)
    params:
        " ".join([get_params("trim_cutadapt_se", "options"),
                  get_params("trim_cutadapt_se", "adapters")])
    threads: get_params("trim_cutadapt_se", "threads")
    log: "logs/{interim}/mapping/trim/{prefix}{fastq}{gz}.log"
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/bio/cutadapt/se"
