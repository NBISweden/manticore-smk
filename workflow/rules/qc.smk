rule qc_all:
    input: multiqc = __REPORTS__ / "qc/multiqc.html",
           fastqc = fastqc_all,
           jellyfish = jellyfish_all

rule multiqc:
    output: "{reports}/qc/multiqc.html"
    input: unpack(multiqc_all)
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['multiqc']['runtime']
    params: config["qc"]["multiqc"]["options"]
    log: "logs/{reports}/qc/multiqc.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/multiqc"


rule fastqc:
    output: html = "{results}/qc/fastqc/{prefix}{bamfastq}{gz}.html",
            zip = "{results}/qc/fastqc/{prefix}{bamfastq}{gz}_fastqc.zip"
    input:  "{prefix}{bamfastq}{gz}"
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['fastqc']['runtime']
    params: config["qc"]["fastqc"]["options"]
    threads: config["qc"]["fastqc"]["threads"]
    log: "logs/{results}/qc/fastqc/{prefix}{bamfastq}{gz}.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/fastqc"


rule jellyfish_count:
    output: jf = temp("{results}/qc/jellyfish/{sample}{trimmed}{ext}{gz}.{kmer}mer_counts.jf")
    input: unpack(jellyfish_count)
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['jellyfish']['count']['runtime']
    wildcard_constraints:
        ext = "(.fq|.fastq|.fa|.fasta|.txt|)",
        trimmed = "(|.trimmed)"
    params:
        options = config['qc']['jellyfish']['count']['options']
    threads: config['qc']['jellyfish']['count']['threads']
    log: "logs/{results}/qc/jellyfish/{sample}{trimmed}{ext}{gz}.{kmer}mer_counts.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/count"


rule jellyfish_histo:
    output: hist = "{prefix}.{kmer}_jf.hist"
    input: counts = "{prefix}.{kmer}mer_counts.jf"
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['jellyfish']['histo']['runtime']
    params:
        options = config['qc']['jellyfish']['histo']['options']
    threads: config['qc']['jellyfish']['histo']['threads']
    log: "logs/{prefix}.{kmer}mer_counts.jf.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/histo"


rule picard_collect_alignment_summary_metrics:
    output: "{results}/qc/align/{prefix}{bam}.align_metrics.txt"
    input: bam = "{prefix}{bam}",
           ref = config['db']['ref']
    log: "logs/{results}/qc/align/{prefix}{bam}.align_metrics.log"
    params:
        config['qc']['picard']['options']
    threads: config['qc']['picard']['collect_alignment_summary_metrics']['threads']
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/collectalignmentsummarymetrics"


rule picard_collect_insert_size_metrics:
    output: txt = "{results}/qc/align/{prefix}{bam}.insert_metrics.txt",
            pdf = "{results}/qc/align/{prefix}{bam}.insert_metrics.pdf",
    input: "{prefix}{bam}"
    log: "logs/{results}/qc/align/{prefix}{bam}.insert_metrics.log"
    params:
        config['qc']['picard']['options']
    threads: config['qc']['picard']['collect_insert_size_metrics']['threads']
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/collectinsertsizemetrics"


rule picard_mark_duplicates:
    output: metrics = "{interim}/map/bwa/{genome}/dedup/{sample}{bam}.dup_metrics.txt",
            bam = "{interim}/map/bwa/{genome}/dedup/{sample}{bam}"
    input: "{interim}/map/bwa/{genome}/{sample}{bam}"
    log: "logs/{interim}/qc/align/{genome}/{sample}{bam}.dup_metrics.log"
    params:
        config['qc']['picard']['mark_duplicates']['options']
    threads: config['qc']['picard']['mark_duplicates']['threads']
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/markduplicates"
