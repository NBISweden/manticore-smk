rule filter_vcf_select:
    """Run filter select step on vcf file"""
    output:
        vcf=temp(
            "{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz"
        ),
        tbi=temp(
            "{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.tbi"
        ),
        cmd="{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.sh",
    input:
        unpack(filter_vcf_input),
        ref=cfg.db.ref,
    wildcard_constraints:
        filtername="select",
    params:
        options=lambda wildcards: cfg.params(wildcards, "options", "filter_vcf_select"),
        java_opts=cfg.ruleconf("filter_vcf_select").java_opts,
    resources:
        runtime=cfg.ruleconf("filter_vcf_select").runtime,
        mem_mb=cfg.ruleconf("filter_vcf_select").mem_mb,
    threads: 1
    log:
        "logs/{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.log",
    wrapper:
        f"{WRAPPER_PREFIX}/bio/filter/vcf_select"


rule filter_vcf_filter:
    """Run filter step on vcf file"""
    output:
        vcf=temp(
            "{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz"
        ),
        tbi=temp(
            "{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.tbi"
        ),
        cmd="{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.sh",
    input:
        unpack(filter_vcf_input),
        ref=cfg.db.ref,
    wildcard_constraints:
        filtername="filter",
    params:
        filters=lambda wildcards: cfg.params(wildcards, "filters"),
        options=lambda wildcards: cfg.params(wildcards, "options", "filter_vcf_filter"),
        java_opts=cfg.ruleconf("filter_vcf_filter").java_opts,
    resources:
        runtime=cfg.ruleconf("filter_vcf_filter").runtime,
        mem_mb=cfg.ruleconf("filter_vcf_filter").mem_mb,
    threads: cfg.ruleconf("filter_vcf_filter").threads
    log:
        "logs/{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.log",
    wrapper:
        f"{WRAPPER_PREFIX}/bio/filter/vcf_filter"


rule filter_vcf_concat:
    """Concatenate/merge vcf files"""
    output:
        vcf="{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz",
        tbi="{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.tbi",
        cmd="{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.sh",
        stats="{results}/qc/variants/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.stats",
    input:
        unpack(filter_vcf_input),
        ref=cfg.db["ref"],
    wildcard_constraints:
        filtername="concat",
    params:
        options=lambda wildcards: cfg.params(wildcards, "options", "filter_vcf_concat"),
        java_opts=cfg.ruleconf("filter_vcf_concat").java_opts,
    resources:
        runtime=cfg.ruleconf("filter_vcf_concat").runtime,
        mem_mb=cfg.ruleconf("filter_vcf_concat").mem_mb,
    threads: cfg.ruleconf("filter_vcf_concat").threads
    log:
        "logs/{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.log",
    wrapper:
        f"{WRAPPER_PREFIX}/bio/filter/vcf_concat"
