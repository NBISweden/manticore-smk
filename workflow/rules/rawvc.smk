rule all_rawvc:
    input: unpack(all_rawvc_input)


rule rawvc_pybedtools_make_bed_targets:
    output:
        bed = "{prefix}/{ind_vc}/{region}.{partition}.bed"
    input: bed = "{prefix}/{region}.bed"
    params:
        npart = lambda wildcards: cfg['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    log: "logs/{prefix}/{ind_vc}/{region}.{partition}.log"
    threads: 1
    wrapper: f"{WRAPPER_PREFIX}/bio/pybedtools/make_bed_targets"


rule rawvc_samtools_faidx_ref:
    output: f"{cfg['db']['ref']}.fai"
    input: cfg['db']['ref']
    log: "logs/db/ref/samtools_faidx.log"
    params:
        cfg.rule("rawvc_samtools_faidx_ref").params("options")
    resources: runtime = lambda wildcards, attempt: cfg.rule("rawvc_samtools_faidx_ref", attempt).resources("runtime")
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/samtools/faidx"


rule rawvc_gatkhc_targets:
    """Run GATK HaplotypeCaller on target file"""
    output:
        vcf = temp("{interim}/{group}/rawvc/gatkhc/{sample}.{target}.{region}{mode}.vcf{gz}"),
        tbi = temp("{interim}/{group}/rawvc/gatkhc/{sample}.{target}.{region}{mode}.vcf{gz}.tbi")
    input:
        unpack(rawvc_gatkhc_targets_input)
    wildcard_constraints: mode = "(.g|)"
    resources:
        runtime = lambda wildcards, attempt: cfg.rule("rawvc_gatkhc_targets", attempt).resources("runtime"),
        mem_mb = lambda wildcards, attempt: cfg.rule("rawvc_gatkhc_targets", attempt).resources( "mem_mb")
    params:
        options = cfg.rule("rawvc_gatkhc_targets").params("options"),
        annotation = cfg.rule("rawvc_gatkhc_targets").params("annotation"),
        ploidy = lambda wildcards: cfg.ploidy(wildcards.sample, wildcards.region)
    threads: cfg.rule("rawvc_gatkhc_targets").threads
    log: "logs/{interim}/{group}/rawvc/gatkhc/{sample}.{target}.{region}{mode}.vcf{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/gatk/hc_targets"


rule rawvc_picard_create_sequence_dictionary:
    """Create sequence dictionary"""
    output: re.sub(wc["fa"], ".dict", cfg['db']['ref'])
    input: cfg['db']['ref']
    params:
        extra = cfg.rule("rawvc_picard_create_sequence_dictionary").params("options")
    threads: cfg.rule("rawvc_picard_create_sequence_dictionary").threads
    log: f"logs/rawvc/picard/{cfg['db']['ref']}.dict"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/createsequencedictionary"


rule rawvc_gatk_genomics_db_import:
    """Import gvcf files to genomics db. By default generates total db and population-based db"""
    output: directory("{results}/{group}/rawvc/gatkhc/genomicsdb/{population}{dot}{region}.{target}.db")
    input: unpack(rawvc_gatk_genomics_db_import_input)
    params:
        options = cfg.rule("rawvc_gatk_genomics_db_import").params("options")
    resources:
        runtime = lambda wildcards, attempt: cfg.rule("rawvc_gatk_genomics_db_import", attempt).resources("runtime"),
        mem_mb = lambda wildcards, attempt: cfg.rule("rawvc_gatk_genomics_db_import", attempt).resources("mem_mb")
    threads: cfg.rule("rawvc_gatk_genomics_db_import").threads
    log: "logs/{results}/{group}/rawvc/gatkhc/genomicsdb/{population}{dot}{region}.{target}.db"
    wrapper: f"{WRAPPER_PREFIX}/bio/gatk/genomics_db_import"


rule rawvc_gatk_genotype_gvcfs:
    output:
        vcf = "{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.{target}.vcf.gz",
        tbi = "{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.{target}.vcf.gz.tbi"
    input: unpack(rawvc_gatk_genotype_gvcfs_input)
    params:
        options = cfg.rule("rawvc_gatk_genotype_gvcfs").params("options"),
        annotation = cfg.rule("rawvc_gatk_genotype_gvcfs").params("annotation")
    resources:
        runtime = lambda wildcards, attempt: cfg.rule("rawvc_gatk_genotype_gvcfs", attempt).resources("runtime"),
        mem_mb = lambda wildcards, attempt: cfg.rule("rawvc_gatk_genotype_gvcfs", attempt).resources("mem_mb")
    threads: cfg.rule("rawvc_gatk_genotype_gvcfs").threads
    log: "logs/{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.{target}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/gatk/genotype_gvcfs"


rule rawvc_bcftools_concat_vcfs_targets:
    """Concatenate target vcf results with bcftools"""
    output:
        vcf = "{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.vcf.gz",
        tbi = "{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.vcf.gz.tbi",
        stats = "{results}/qc/variants/{group}/rawvc/gatkhc/{population}{dot}{region}.vcf.gz.stats"
    input: unpack(rawvc_bcftools_concat_vcfs_targets_input),
           ref = cfg.db.ref
    params:
        extra = cfg.rule("rawvc_bcftools_concat_vcfs_targets").params("options")
    resources:
        runtime = lambda wildcards, attempt: cfg.rule("rawvc_bcftools_concat_vcfs_targets", attempt).resources("runtime"),
        mem_mb = lambda wildcards, attempt: cfg.rule("rawvc_bcftools_concat_vcfs_targets", attempt).resources("mem_mb")
    threads: cfg.rule("rawvc_bcftools_concat_vcfs_targets").threads
    log: "logs/{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/bcftools/concat"


localrules: rawvc_pybedtools_make_bed_targets, rawvc_picard_create_sequence_dictionary
