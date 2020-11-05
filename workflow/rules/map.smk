rule all_map:
    input: all_bwa_mem_samples

rule bwa_link_ref:
    output: "{interim}/map/bwa/index/{genome}.fasta"
    input: ref = config["db"]["ref"]
    params:
        seq = lambda wildcards, input: os.path.abspath(input.ref)
    log: "logs/{interim}/map/bwa/index/{genome}.fasta.log"
    shell: "ln -s {params.seq} {output}"


rule map_samtools_faidx_ref:
    output: "{interim}/map/{prefix}{fa}{bgz}.fai"
    input: "{interim}/map/{prefix}{fa}{bgz}"
    log: "logs/{interim}/map/{prefix}{fa}{bgz}.log"
    params: config['map']['samtools']['faidx']['options']
    resources: runtime = lambda wilcards, attempt: attempt * config['map']['samtools']['faidx']['runtime']
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/samtools/faidx"


rule map_samtools_index:
    output: "{interim}/map/{prefix}{bam}.bai"
    input: "{interim}/map/{prefix}{bam}"
    log: "logs/{interim}/map/{prefix}{bam}.log"
    params: config['map']['samtools']['index']['options']
    resources: runtime = lambda wilcards, attempt: attempt * config['map']['samtools']['index']['runtime']
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/samtools/index"


rule bwa_index:
    """bwa index a reference"""
    output: expand("{{interim}}/map/bwa/index/{{genome}}{ext}", ext=bwa_ext)
    input: "{interim}/map/bwa/index/{genome}.fasta"
    resources: runtime = lambda wildcards, attempt: attempt * config['map']['bwa']['index']['runtime']
    log: "logs/{interim}/map/bwa/index/{genome}.log"
    params:
        prefix = lambda wildcards: f"{wildcards.interim}/map/bwa/index/{wildcards.genome}",
        algorithm = "bwtsw"
    threads: config["map"]["bwa"]["index"]["threads"]
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/bwa/index"


rule bwa_mem:
    output: temp("{interim}/map/bwa/{sample}/{unit}.bam")
    input:
        index = bwa_mem_index_ext,
        reads = bwa_mem_input
    resources: runtime = lambda wildcards, attempt: attempt * config['map']['bwa']['mem']['runtime']
    log: "logs/{interim}/map/bwa/{sample}/{unit}.log"
    params:
        index = bwa_mem_index,
        sort = config["map"]["bwa"]["mem"]["sort"],
        sort_order = config["map"]["bwa"]["mem"]["sort_order"],
        sort_extra = config["map"]["bwa"]["mem"]["sort_extra"],
        extra = bwa_mem_rg
    threads: config["map"]["bwa"]["mem"]["threads"]
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/bwa/mem"


rule picard_merge_sam:
    output: "{interim}/map/bwa/{sample}.bam"
    input: picard_merge_sam_input
    resources: runtime = lambda wildcards, attempt: attempt * config['map']['bwa']['mem']['runtime']
    log: "logs/{interim}/map/bwa/{sample}.log"
    threads: config["map"]["bwa"]["mem"]["threads"]
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/mergesamfiles"

localrules: bwa_link_ref
