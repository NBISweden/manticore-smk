rule pybedtools_make_bed_targets:
    output:
        bed = temp("{prefix}/{ind_vc}/{region}.{partition}.bed")
    input: bed = "{prefix}/{region}.bed"
    params:
        npart = lambda wildcards: config['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    log: "logs/{prefix}/{ind_vc}/{region}.{partition}.log"
    threads: 1
    wrapper: f"{WRAPPER_PREFIX}/bio/pybedtools/make_bed_targets"


rule rawvc_samtools_faidx_ref:
    output: f"{config['db']['ref']}.fai"
    input: config['db']['ref']
    log: "logs/db/ref/samtools_faidx.log"
    params: config['rawvc']['samtools']['faidx']['options']
    resources: runtime = lambda wilcards, attempt: attempt * config['rawvc']['samtools']['faidx']['runtime']
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/samtools/faidx"


rule gatkhc_targets:
    """Run GATK HaplotypeCaller on target file"""
    output:
        vcf = temp("{interim}/rawvc/gatkhc/{sample}.{target}.{region}{mode}.vcf{bgz}"),
        tbi = temp("{interim}/rawvc/gatkhc/{sample}.{target}.{region}{mode}.vcf{bgz}.tbi")
    input:
        unpack(gatkhc_targets_input)
    wildcard_constraints: mode = "(.g|)"
    params: options = config['rawvc']['gatk']['haplotype_caller_targets']['options']
    threads: config['rawvc']['gatk']['haplotype_caller_targets']['threads']
    log: "logs/{interim}/rawvc/gatkhc/{sample}.{target}.{region}{mode}.vcf{bgz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/gatk/hc_targets"


rule rawvc_picard_create_sequence_dictionary:
    """Create sequence dictionary"""
    output: re.sub(wc["fa"], ".dict", config['db']['ref'])
    input: config['db']['ref']
    params: extra = config['rawvc']['picard']['create_sequence_dictionary']['options']
    threads: config['rawvc']['picard']['create_sequence_dictionary']['threads']
    log: f"logs/rawvc/picard/{config['db']['ref']}.dict"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/createsequencedictionary"



localrules: pybedtools_make_bed_targets
