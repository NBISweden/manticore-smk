rule all_map:
    """Target to finish all mapping, including merging and deduplication"""
    input: all_map_input

rule map_bwa_link_ref:
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
    params: get_params("map_samtools_faidx_ref", "options")
    resources:
        runtime = lambda wilcards, attempt: resources("map_samtools_faidx_ref", "runtime", attempt)
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/samtools/faidx"


rule map_samtools_index:
    output: "{interim}/map/{prefix}{bam}.bai"
    input: "{interim}/map/{prefix}{bam}"
    log: "logs/{interim}/map/{prefix}{bam}.log"
    params: get_params("map_samtools_index", "options")
    resources:
        runtime = lambda wilcards, attempt: resources("map_samtools_index", "runtime", attempt)
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/samtools/index"


rule map_bwa_index:
    """bwa index a reference"""
    output: expand("{{interim}}/map/bwa/index/{{genome}}{ext}", ext=bwa_ext)
    input: "{interim}/map/bwa/index/{genome}.fasta"
    resources:
        runtime = lambda wilcards, attempt: resources("map_bwa_index", "runtime", attempt)
    log: "logs/{interim}/map/bwa/index/{genome}.log"
    params:
        prefix = lambda wildcards: f"{wildcards.interim}/map/bwa/index/{wildcards.genome}",
        algorithm = "bwtsw"
    threads: get_params("map_bwa_index", "threads")
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/bwa/index"


rule map_bwa_mem:
    output: temp("{interim}/map/bwa/{sample}/{unit}.bam")
    input:
        index = bwa_mem_index_ext,
        reads = bwa_mem_input
    resources:
        runtime = lambda wilcards, attempt: resources("map_bwa_mem", "runtime", attempt)
    log: "logs/{interim}/map/bwa/{sample}/{unit}.log"
    params:
        index = bwa_mem_index,
        sort = get_params("map_bwa_mem", "sort"),
        sort_order = get_params("map_bwa_mem", "sort_order"),
        sort_extra = get_params("map_bwa_mem", "sort_extra"),
        extra = bwa_mem_rg
    threads: get_params("map_bwa_mem", "threads")
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/bwa/mem"


rule map_picard_merge_sam:
    output: "{interim}/map/bwa/{sample}.bam"
    input: picard_merge_sam_input
    resources:
        runtime = lambda wildcards, attempt: resources("map_picard_merge_sam", "runtime", attempt)
    log: "logs/{interim}/map/bwa/{sample}.log"
    threads: get_params("map_picard_merge_sam", "threads")
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/mergesamfiles"

localrules: map_bwa_link_ref
