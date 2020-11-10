# Rules for manual variant filtering

rule all_variantfilter:
    input: []


rule filter_gatk_select_variants:
    """Run GATK SelectVariants to select SNPs"""
    output: vcf = "{results}/{callset}/{caller}/select/{region}.{vartype}.bcf",
            idx = "{results}/{callset}/{caller}/select/{region}.{vartype}.bcf.idx"
    input: vcf = "{results}/{callset}/{caller}/unfiltered/{region}.bcf",
           idx = "{results}/{callset}/{caller}/unfiltered/{region}.bcf.idx",
           ref = config['db']['ref']
    wildcard_constraints:
        vartype = "(snp|indel)"
    params:
        extra = lambda wildcards: "--select-type-to-include SNP" if wildcards.vartype == "snp" else "--select-type-to-include INDEL",
        java_opts = get_params("filter_gatk_select_variants", "java_options")
    resources:
        runtime = lambda wildcards, attempt: resources("filter_gatk_select_variants", "runtime"),
        mem_mb = lambda wildcards, attempt: resources("filter_gatk_select_variants", "mem_mb")
    threads: get_params("filter_gatk_select_variants", "threads")
    log: "logs/{results}/{callset}/{caller}/select/{region}.{vartype}.bcf.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/gatk/selectvariants"


rule filter_gatk_variant_JEXL_filtration:
    """Run GATK VariantFiltration using hard filters

    Perform hard filtering using JEXL expressions
    """
    output: vcf = "{results}/{callset}/{caller}/filter/{region}.{vartype}.jexl.bcf"
    input: vcf = "{results}/{callset}/{caller}/select/{region}.{vartype}.bcf",
           idx = "{results}/{callset}/{caller}/select/{region}.{vartype}.bcf.idx",
           ref = config['db']['ref']
    wildcard_constraints:
        vartype = "(snp|indel)"
    params:
        filters = lambda wildcards: get_variant_filters("filter_gatk_variant_JEXL_filtration", wildcards),
        extra = get_params("filter_gatk_variant_JEXL_filtration", "options"),
        java_opts = get_params("filter_gatk_variant_JEXL_filtration", "java_options")
    resources:
        runtime = lambda wildcards, attempt: resources("filter_gatk_variant_JEXL_filtration", "runtime"),
        mem_mb = lambda wildcards, attempt: resources("filter_gatk_variant_JEXL_filtration", "mem_mb")
    threads: get_params("filter_gatk_variant_JEXL_filtration", "threads")
    log: "logs/{results}/{callset}/{caller}/select/{region}.{vartype}.bcf.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/gatk/variantfiltration"
