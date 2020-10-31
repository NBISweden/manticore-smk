rule map_all:
    input: bwa_mem_all

rule bwa_link_ref:
    output: "{interim}/map/bwa/index/{genome}.fasta"
    input: lambda wildcards: config["db"]["ref"]["seq"]
    params:
        seq = lambda wildcards, input: os.path.abspath(input[0])
    shell: "ln -s {params.seq} {output}"

rule bwa_index:
    """bwa index a reference"""
    output: expand("{{interim}}/map/bwa/index/{{genome}}{ext}", ext=bwa_ext)
    input: "{interim}/map/bwa/index/{genome}.fasta"
    resources: runtime = lambda wildcards, attempt: attempt * config['qc']['bwa']['index']['runtime']
    log: "logs/{interim}/map/bwa/index/{genome}.log"
    params:
        prefix = "{interim}/map/bwa/index/{genome}",
        algorithm = "bwtsw"
    threads: config["map"]["bwa"]["index"]["threads"]
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/bwa/index"


rule bwa_mem:
    output: temp("{interim}/map/bwa/{genome}/{sample}/{unit}.bam")
    input:
        index = expand("{{interim}}/map/bwa/index/{{genome}}{ext}", ext=bwa_ext),
        reads = bwa_mem_input
    resources: runtime = lambda wildcards, attempt: attempt * config['qc']['bwa']['mem']['runtime']
    log: "logs/{interim}/map/bwa/{genome}/{sample}/{unit}.log"
    params:
        index = "{interim}/map/bwa/index/{genome}",
        sort = config["map"]["bwa"]["mem"]["sort"],
        sort_order = config["map"]["bwa"]["mem"]["sort_order"],
        sort_extra = config["map"]["bwa"]["mem"]["sort_extra"],
        extra = bwa_mem_rg
    threads: config["map"]["bwa"]["mem"]["threads"]
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/bwa/mem"


rule picard_merge_sam:
    output: "{interim}/map/bwa/{genome}/{sample}.bam"
    input: picard_merge_sam_input
    resources: runtime = lambda wildcards, attempt: attempt * config['map']['bwa']['mem']['runtime']
    log: "logs/{interim}/map/bwa/{genome}/{sample}.log"
    threads: config["map"]["bwa"]["mem"]["threads"]
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/mergesamfiles"

localrules: bwa_link_ref
