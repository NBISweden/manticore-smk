rule filter_vcf_select:
    """Run filter select step on vcf file"""
    output:
        vcf = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz"),
        tbi = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.tbi"),
        cmd = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.sh"
    input:
        unpack(filter_vcf_input),
        ref = cfg.db.ref
    wildcard_constraints:
        filtername = "select"
    params:
        options = lambda wildcards: get_filter_options(wildcards),
        java_opts = lambda wildcards: cfg.rule("filter_vcf_select").params("java_options")
    resources:
        runtime = lambda wildcards, attempt: cfg.rule("filter_vcf_select", attempt).resources("runtime"),
        mem_mb = lambda wildcards, attempt: cfg.rule("filter_vcf_select", attempt).resources("mem_mb"),
    threads: cfg.rule("filter_vcf_select").threads
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/filter/vcf_select"


rule filter_vcf_filter:
    """Run filter step on vcf file"""
    output:
        vcf = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz"),
        tbi = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.tbi"),
        cmd = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.sh"
    input:
        unpack(filter_vcf_input),
        ref = cfg.db['ref']
    wildcard_constraints:
        filtername = "filter"
    params:
        filters = lambda wildcards: get_filter_options(wildcards, key="filters"),
        options = lambda wildcards: get_filter_options(wildcards),
        java_opts = lambda wildcards: cfg.rule("filter_vcf_filter").params("java_options")
    resources:
        runtime = lambda wildcards, attempt: cfg.rule("filter_vcf_filter", attempt).resources("runtime"),
        mem_mb = lambda wildcards, attempt: cfg.rule("filter_vcf_filter", attempt).resources("mem_mb"),
    threads: cfg.rule("filter_vcf_filter").threads
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/filter/vcf_filter"


rule filter_vcf_concat:
    """Concatenate/merge vcf files"""
    output:
        vcf = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz",
        tbi = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.tbi",
        cmd = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.sh",
        stats = "{results}/qc/variants/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.stats"
    input:
        unpack(filter_vcf_input),
        ref = cfg.db['ref']
    wildcard_constraints:
        filtername = "concat"
    params:
        options = lambda wildcards: get_filter_options(wildcards),
        java_opts = lambda wildcards: cfg.rule("filter_vcf_concat").params("java_options")
    resources:
        runtime = lambda wildcards, attempt: cfg.rule("filter_vcf_concat", attempt).resources("runtime"),
        mem_mb = lambda wildcards, attempt: cfg.rule("filter_vcf_concat", attempt).resources("mem_mb"),
    threads: cfg.rule("filter_vcf_concat").threads
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/filter/vcf_concat"

# rule filter_gatk_select_variants:
#     """Run analysis filter"""
#     output:
#         vcf = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz"),
#         tbi = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.tbi")
#     input:
#         unpack(filter_gatk_select_variants_input),
#         ref = cfg.db['ref']
#     wildcard_constraints:
#         filternum = "[0-9]{2}",
#         #tool = "gatk_select_variants"
#     params:
#         extra = lambda wildcards: get_filter_options(wildcards),
#         java_opts = lambda wildcards: cfg.rule("filter_gatk_select_variants").params("java_options")
#     resources:
#         runtime = lambda wildcards, attempt: cfg.rule("filter_gatk_select_variants", attempt).resources("runtime"),
#         mem_mb = lambda wildcards, attempt: cfg.rule("filter_gatk_select_variants", attempt).resources("mem_mb")
#     threads: cfg.rule("filter_gatk_select_variants").threads
#     log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.log"
#     wrapper: f"{SMK_WRAPPER_PREFIX}/bio/gatk/selectvariants"


# rule filter_gatk_jexl_filter_variants:
#     """Run analysis filter"""
#     output:
#         vcf = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz"),
#         tbi = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.tbi")
#     input:
#         unpack(filter_gatk_jexl_filter_variants_input),
#         ref = cfg.db['ref']
#     wildcard_constraints:
#         filternum = "[0-9]{2}",
#         #tool = "gatk_jexl_filter_variants"
#     params:
#         filters = lambda wildcards: get_filter_options(wildcards, gatk=True),
#         extra = lambda wildcards: cfg.rule("filter_gatk_jexl_filter_variants").params("options"),
#         java_opts = lambda wildcards: cfg.rule("filter_gatk_jexl_filter_variants").params("java_options")
#     resources:
#         runtime = lambda wildcards, attempt: cfg.rule("filter_gatk_jexl_filter_variants", attempt).resources("runtime"),
#         mem_mb = lambda wildcards, attempt: cfg.rule("filter_gatk_jexl_filter_variants", attempt).resources("mem_mb")
#     threads: cfg.rule("filter_gatk_select_variants").threads
#     log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.{target}.vcf.gz.log"
#     wrapper: f"{SMK_WRAPPER_PREFIX}/bio/gatk/variantfiltration"


# rule filter_bcftools_concat_vcfs:
#     """Concatenate vcf results with bcftools"""
#     output:
#         vcf = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz",
#         tbi = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.tbi",
#         stats = "{results}/qc/variants/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.stats"
#     input: unpack(filter_bcftools_concat_vcfs_input),
#            ref = cfg.db.ref
#     wildcard_constraints:
#         filternum = "[0-9]{2}",
#         #tool = "bcftools_concat"
#     params:
#         extra = cfg.rule("filter_bcftools_concat_vcfs").params("options")
#     resources:
#         runtime = lambda wildcards, attempt: cfg.rule("filter_bcftools_concat_vcfs", attempt).resources("runtime"),
#         mem_mb = lambda wildcards, attempt: cfg.rule("filter_bcftools_concat_vcfs", attempt).resources("mem_mb")
#     threads: cfg.rule("filter_bcftools_concat_vcfs").threads
#     log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{population}{dot}{region}.vcf.gz.log"
#     wrapper: f"{WRAPPER_PREFIX}/bio/bcftools/concat"
