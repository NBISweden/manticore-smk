rule filter_gatk_select_variants:
    """Run analysis filter"""
    output:
        vcf = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.{target}.vcf.gz"),
        tbi = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.{target}.vcf.gz.tbi")
    input:
        unpack(filter_gatk_select_variants_input),
        ref = config['db']['ref']
    wildcard_constraints:
        filternum = "[0-9]{2}",
        filtertype = "gatk_select_variants"
    params:
        extra = lambda wildcards: get_filter_options(wildcards),
        java_opts = lambda wildcards: get_params("filter_gatk_select_variants", "java_options")
    resources:
        runtime = lambda wildcards, attempt: resources("filter_gatk_select_variants", "runtime"),
        mem_mb = lambda wildcards, attempt: resources("filter_gatk_select_variants", "mem_mb")
    threads: get_params("filter_gatk_select_variants", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.{target}.vcf.gz.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/gatk/selectvariants"


rule filter_gatk_jexl_filter_variants:
    """Run analysis filter"""
    output:
        vcf = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.{target}.vcf.gz"),
        tbi = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.{target}.vcf.gz.tbi")
    input:
        unpack(filter_gatk_jexl_filter_variants_input),
        ref = config['db']['ref']
    wildcard_constraints:
        filternum = "[0-9]{2}",
        filtertype = "gatk_jexl_filter_variants"
    params:
        filters = lambda wildcards: get_filter_options(wildcards, gatk=True),
        extra = lambda wildcards: get_params("filter_gatk_jexl_filter_variants", "options"),
        java_opts = lambda wildcards: get_params("filter_gatk_jexl_filter_variants", "java_options")
    resources:
        runtime = lambda wildcards, attempt: resources("filter_gatk_jexl_filter_variants", "runtime"),
        mem_mb = lambda wildcards, attempt: resources("filter_gatk_jexl_filter_variants", "mem_mb")
    threads: get_params("filter_gatk_select_variants", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.{target}.vcf.gz.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/gatk/variantfiltration"


rule filter_bcftools_concat_vcfs:
    """Concatenate vcf results with bcftools"""
    output:
        vcf = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.vcf.gz",
        tbi = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.vcf.gz.tbi",
        stats = "{results}/qc/variants/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.vcf.gz.stats"
    input: unpack(filter_bcftools_concat_vcfs_input),
           ref = config["db"]["ref"]
    wildcard_constraints:
        filternum = "[0-9]{2}",
        filtertype = "bcftools_concat"
    params:
        extra = get_params("filter_bcftools_concat_vcfs", "options")
    resources:
        runtime = lambda wildcards, attempt: resources("filter_bcftools_concat_vcfs", "runtime"),
        mem_mb = lambda wildcards, attempt: resources("filter_bcftools_concat_vcfs", "mem_mb")
    threads: get_params("filter_bcftools_concat_vcfs", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/bcftools/concat"
