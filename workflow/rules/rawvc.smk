rule all_rawvc:
    input:
        unpack(all_rawvc_input),


rule rawvc_pybedtools_make_bed_targets:
    output:
        bed="{prefix}/{ind_vc}/{region}.{partition}.bed",
    input:
        bed="{prefix}/{region}.bed",
    params:
        npart=lambda wildcards: cfg["workflow"]["regions"]
        .get(wildcards.region, {})
        .get("npart", 1),
    log:
        "logs/{prefix}/{ind_vc}/{region}.{partition}.log",
    threads: 1
    wrapper:
        f"{WRAPPER_PREFIX}/bio/pybedtools/make_bed_targets"


rule rawvc_samtools_faidx_ref:
    output:
        f"{cfg['db']['ref']}.fai",
    input:
        cfg["db"]["ref"],
    log:
        "logs/db/ref/samtools_faidx.log",
    params:
        cfg.ruleconf("rawvc_samtools_faidx_ref").options,
    resources:
        runtime=cfg.ruleconf("rawvc_samtools_faidx_ref").runtime,
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/bio/samtools/faidx"


rule rawvc_gatkhc_targets:
    """Run GATK HaplotypeCaller on target file"""
    output:
        vcf=temp(
            "{interim}/{group}/rawvc/gatkhc/{sample}.{target}.{region}{mode}.vcf{gz}"
        ),
        tbi=temp(
            "{interim}/{group}/rawvc/gatkhc/{sample}.{target}.{region}{mode}.vcf{gz}.tbi"
        ),
    input:
        unpack(rawvc_gatkhc_targets_input),
    wildcard_constraints:
        mode="(.g|)",
    resources:
        runtime=cfg.ruleconf("rawvc_gatkhc_targets").runtime,
        mem_mb=cfg.ruleconf("rawvc_gatkhc_targets").mem_mb,
    params:
        options=cfg.ruleconf("rawvc_gatkhc_targets").options,
        annotation=cfg.ruleconf("rawvc_gatkhc_targets").extra["annotation"],
        ploidy=lambda wildcards: cfg.ploidy(wildcards.region, sample=wildcards.sample),
    threads: cfg.ruleconf("rawvc_gatkhc_targets").threads
    log:
        "logs/{interim}/{group}/rawvc/gatkhc/{sample}.{target}.{region}{mode}.vcf{gz}.log",
    wrapper:
        f"{WRAPPER_PREFIX}/bio/gatk/hc_targets"


rule rawvc_picard_create_sequence_dictionary:
    """Create sequence dictionary"""
    output:
        re.sub(wc["fa"], ".dict", cfg.db.ref),
    input:
        cfg.db.ref,
    params:
        extra=cfg.ruleconf("rawvc_picard_create_sequence_dictionary").options,
    threads: cfg.ruleconf("rawvc_picard_create_sequence_dictionary").threads
    log:
        f"logs/rawvc/picard/{cfg['db']['ref']}.dict.log",
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/bio/picard/createsequencedictionary"


rule rawvc_gatk_genomics_db_import:
    """Import gvcf files to genomics db. By default generates total db and population-based db"""
    output:
        directory(
            "{results}/{group}/rawvc/gatkhc/genomicsdb/{population}{dot}{region}.{target}.db"
        ),
    input:
        unpack(rawvc_gatk_genomics_db_import_input),
    params:
        options=cfg.ruleconf("rawvc_gatk_genomics_db_import").options,
    resources:
        runtime=cfg.ruleconf("rawvc_gatk_genomics_db_import").runtime,
        mem_mb=cfg.ruleconf("rawvc_gatk_genomics_db_import").mem_mb,
    threads: cfg.ruleconf("rawvc_gatk_genomics_db_import").threads
    log:
        "logs/{results}/{group}/rawvc/gatkhc/genomicsdb/{population}{dot}{region}.{target}.db",
    wrapper:
        f"{WRAPPER_PREFIX}/bio/gatk/genomics_db_import"


rule rawvc_gatk_genotype_gvcfs:
    output:
        vcf="{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.{target}.vcf.gz",
        tbi="{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.{target}.vcf.gz.tbi",
    input:
        unpack(rawvc_gatk_genotype_gvcfs_input),
    params:
        options=cfg.ruleconf("rawvc_gatk_genotype_gvcfs").options,
        annotation=cfg.ruleconf("rawvc_gatk_genotype_gvcfs").extra["annotation"],
    resources:
        runtime=cfg.ruleconf("rawvc_gatk_genotype_gvcfs").runtime,
        mem_mb=cfg.ruleconf("rawvc_gatk_genotype_gvcfs").mem_mb,
    threads: cfg.ruleconf("rawvc_gatk_genotype_gvcfs").threads
    log:
        "logs/{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.{target}.vcf.gz.log",
    wrapper:
        f"{WRAPPER_PREFIX}/bio/gatk/genotype_gvcfs"


rule rawvc_bcftools_concat_vcfs_targets:
    """Concatenate target vcf results with bcftools"""
    output:
        vcf="{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.vcf.gz",
        tbi="{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.vcf.gz.tbi",
        stats="{results}/qc/variants/{group}/rawvc/gatkhc/{population}{dot}{region}.vcf.gz.stats",
    input:
        unpack(rawvc_bcftools_concat_vcfs_targets_input),
        ref=cfg.db.ref,
    params:
        extra=cfg.ruleconf("rawvc_bcftools_concat_vcfs_targets").options,
    resources:
        runtime=cfg.ruleconf("rawvc_bcftools_concat_vcfs_targets").runtime,
        mem_mb=cfg.ruleconf("rawvc_bcftools_concat_vcfs_targets").mem_mb,
    threads: cfg.ruleconf("rawvc_bcftools_concat_vcfs_targets").threads
    log:
        "logs/{results}/{group}/rawvc/gatkhc/{population}{dot}{region}.vcf.gz.log",
    wrapper:
        f"{WRAPPER_PREFIX}/bio/bcftools/concat"


localrules:
    rawvc_pybedtools_make_bed_targets,
    rawvc_picard_create_sequence_dictionary,
