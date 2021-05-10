rule filter_vcf_select:
    """Run filter select step on vcf file"""
    output:
        vcf = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.{target}.vcf.gz"),
        tbi = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.{target}.vcf.gz.tbi"),
        cmd = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.{target}.vcf.gz.sh"
    input:
        unpack(filter_vcf_input),
        ref = config['db']['ref']
    wildcard_constraints:
        filtername = "select"
    params:
        options = lambda wildcards: get_filter_options(wildcards),
        java_opts = lambda wildcards: get_params("filter_vcf_select", "java_options")
    resources:
        runtime = lambda wildcards, attempt: resources("filter_vcf_select", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("filter_vcf_select", "mem_mb", attempt),
    threads: get_params("filter_vcf_select", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.{target}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/filter/vcf_select"


rule filter_vcf_filter:
    """Run filter step on vcf file"""
    output:
        vcf = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.{target}.vcf.gz"),
        tbi = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.{target}.vcf.gz.tbi"),
        cmd = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.{target}.vcf.gz.sh"
    input:
        unpack(filter_vcf_input),
        ref = config['db']['ref']
    wildcard_constraints:
        filtername = "filter"
    params:
        filters = lambda wildcards: get_filter_options(wildcards, key="filters"),
        options = lambda wildcards: get_filter_options(wildcards),
        java_opts = lambda wildcards: get_params("filter_vcf_filter", "java_options")
    resources:
        runtime = lambda wildcards, attempt: resources("filter_vcf_filter", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("filter_vcf_filter", "mem_mb", attempt),
    threads: get_params("filter_vcf_filter", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.{target}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/filter/vcf_filter"


rule filter_vcf_concat:
    """Concatenate/merge vcf files"""
    output:
        vcf = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.vcf.gz",
        tbi = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.vcf.gz.tbi",
        cmd = "{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.vcf.gz.sh",
        stats = "{results}/qc/variants/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.vcf.gz.stats"
    input:
        unpack(filter_vcf_input),
        ref = config['db']['ref']
    wildcard_constraints:
        filtername = "concat"
    params:
        options = lambda wildcards: get_filter_options(wildcards),
        java_opts = lambda wildcards: get_params("filter_vcf_concat", "java_options")
    resources:
        runtime = lambda wildcards, attempt: resources("filter_vcf_concat", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("filter_vcf_concat", "mem_mb", attempt),
    threads: get_params("filter_vcf_concat", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{tool}/{region}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/filter/vcf_concat"

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
        runtime = lambda wildcards, attempt: resources("filter_gatk_select_variants", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("filter_gatk_select_variants", "mem_mb", attempt)
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
        runtime = lambda wildcards, attempt: resources("filter_gatk_jexl_filter_variants", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("filter_gatk_jexl_filter_variants", "mem_mb", attempt)
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
        runtime = lambda wildcards, attempt: resources("filter_bcftools_concat_vcfs", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("filter_bcftools_concat_vcfs", "mem_mb", attempt)
    threads: get_params("filter_bcftools_concat_vcfs", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{region}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/bcftools/concat"
