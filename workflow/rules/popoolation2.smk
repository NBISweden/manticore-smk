rule all_popoolation2:
    """Run popoolation2 analyses"""
    input: unpack(all_popoolation2_input)


rule all_popoolation2_raw:
    """Generate raw popoolation2 data"""
    input: unpack(all_popoolation2_raw_input)


# Rules
rule popoolation2_samtools_mpileup:
    """Run samtools mpileup for popoolation2"""
    output: pileup = temp("{interim_pool}/raw/popoolation2/{sex}.{region}.{target}.mpileup{gz}")
    input: unpack(popoolation2_samtools_mpileup_input)
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation2_samtools_mpileup", "runtime", attempt)
    params:
        options = get_params("popoolation2_samtools_mpileup", "options")
    threads: lambda wildcards: resources("popoolation2_samtools_mpileup", "threads")
    log: "logs/{interim_pool}/raw/popoolation2/{sex}.{region}.{target}.mpileup{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/samtools_mpileup"


rule popoolation2_mpileup2sync_jar:
    """Convert mpileup output to sync format"""
    output: sync = "{results}/{group}/raw/popoolation2/{sex}.{region}.{target}.sync{gz}"
    input: mpileup = __INTERIM__ / "{group}/raw/popoolation2/{sex}.{region}.{target}.mpileup.gz"
    resources:
        mem_mb = lambda wildcards, attempt: resources("popoolation2_mpileup2sync_jar", "mem_mb", attempt),
        runtime = lambda wildcards, attempt: resources("popoolation2_mpileup2sync_jar", "runtime", attempt)
    params:
        java_options = get_params("popoolation2_mpileup2sync_jar", "java_options")
    threads: get_params("popoolation2_mpileup2sync_jar", "threads")
    log: "logs/{results}/{group}/raw/popoolation2/{sex}.{region}.{target}.sync{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/mpileup2sync_jar"


rule popoolation2_indel_filtering_identify_indel_regions:
    """Identify indel regions in mpileup file"""
    output: gtf = "{results_pool}/raw/popoolation2.indels/{sex}.{region}.{target}.mpileup{gz}.indels.gtf"
    input: mpileup = __INTERIM_POOL__ / "raw/popoolation2/{sex}.{region}.{target}.mpileup{gz}"
    resources:
        mem_mb = lambda wildcards, attempt: resources("popoolation2_indel_filtering_identify_indel_regions", "mem_mb", attempt),
        runtime = lambda wildcards, attempt: resources("popoolation2_indel_filtering_identify_indel_regions", "runtime", attempt)
    params:
        java_options = get_params("popoolation2_indel_filtering_identify_indel_regions", "java_options")
    threads: get_params("popoolation2_indel_filtering_identify_indel_regions", "threads")
    log: "logs/{results_pool}/raw/popoolation2.indels/{sex}.{region}.{target}.mpileup{gz}.indels.gtf.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/indel_filtering_identify_indel_regions"


rule popoolation2_indel_filtering_filter_sync_by_gtf:
    """Filter sync file by indel gtf input generated by popoolation2_indel_filtering_identify_indel_regions"""
    output: sync = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sex}.{region}.{target}.sync{gz}")
    input: unpack(get_popoolation2_filter_input),
           #sync = "{results}/{group}/raw/popoolation2/{sex}.{region}.{target}.sync{gz}",
           gtf = popoolation2_filter_pileup_by_gtf_input
    wildcard_constraints:
        filternum = "[0-9]{2}",
        filtertype = "popoolation2_indel_filtering_filter_sync_by_gtf"
    resources:
        mem_mb = lambda wildcards, attempt: resources("popoolation2_indel_filtering_filter_sync_by_gtf", "mem_mb", attempt),
        runtime = lambda wildcards, attempt: resources("popoolation2_indel_filtering_filter_sync_by_gtf", "runtime", attempt)
    params:
        java_options = get_params("popoolation2_indel_filtering_filter_sync_by_gtf", "java_options")
    threads: get_params("popoolation2_indel_filtering_filter_sync_by_gtf", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sex}.{region}.{target}.sync{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/indel_filtering_filter_sync_by_gtf"


rule popoolation2_indel_filtering_remove_indels:
    """Remove positions with indels"""
    output: sync = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sex}.{region}.{target}.sync{gz}")
    input: unpack(get_popoolation2_filter_input)
    wildcard_constraints:
        filternum = "[0-9]{2}",
        filtertype = "popoolation2_indel_filtering_remove_indels"
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation2_indel_filtering_remove_indels", "runtime", attempt)
    params:
        options = get_params("popoolation2_indel_filtering_remove_indels", "options")
    threads: get_params("popoolation2_indel_filtering_remove_indels", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sex}.{region}.{target}.sync{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/indel_filtering_remove_indels"


rule popoolation2_subsample_synchronized:
    """Subsample synchronized reads"""
    output:
        sync = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sex}.{region}.{target}.sync{gz}")
    input: unpack(get_popoolation2_filter_input)
    wildcard_constraints:
        filternum = "[0-9]{2}",
        filtertype = "popoolation2_subsample_synchronized"
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation2_subsample_synchronized", "runtime", attempt)
    params:
        options = lambda wildcards: get_filter_options(wildcards)
    threads: get_params("popoolation2_subsample_synchronized", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sex}.{region}.{target}.sync{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/subsample_synchronized"


rule popoolation2_fst_sliding:
    """Calculate Fst values using a sliding window approach"""
    output: fst = "{results}/{group}/analysis/{analysis}/{statname}/{sex}.{region}.w{window_size}.s{step_size}.{target}.fst.txt.gz"
    input: unpack(get_popoolation2_filter_input)
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation2_fst_sliding", "runtime", attempt)
    params:
        options = lambda wildcards: get_stat_options(wildcards, rulename="popoolation2_fst_sliding"),
        samples = lambda wildcards: pools if wildcards.sex == "common" else pools[pools.sex.isin([wildcards.sex])]
    threads: lambda wildcards: resources("popoolation2_fst_sliding", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{statname}/{sex}.{region}.w{window_size}.s{step_size}.{target}.fst.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/fst_sliding"


rule popoolation2_fisher_test:
    """Run Fisher exact test to estimate the significance of allele frequency differences"""
    output: fet = "{results}/{group}/analysis/{analysis}/{statname}/{sex}.{region}.w{window_size}.s{step_size}.{target}.fet.txt.gz"
    input: unpack(get_popoolation2_filter_input)
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation2_fisher_test", "runtime", attempt)
    params:
        options = lambda wildcards: get_stat_options(wildcards, rulename="popoolation2_fisher_test"),
        samples = lambda wildcards: pools if wildcards.sex == "common" else pools[pools.sex.isin([wildcards.sex])]
    threads: lambda wildcards: resources("popoolation2_fisher_test", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{statname}/{sex}.{region}.w{window_size}.s{step_size}.{target}.fet.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/fisher_test"


rule popoolation2_snp_frequency_diff:
    """Calculate allele frequency differences"""
    output:
        rc = "{results}/{group}/analysis/{analysis}/{statname}/{sex}.{region}.{target}.sync_rc.txt{gz}",
        pwc = "{results}/{group}/analysis/{analysis}/{statname}/{sex}.{region}.{target}.sync_pwc.txt{gz}"
    input: unpack(get_popoolation2_filter_input)
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation2_snp_frequency_diff", "runtime", attempt)
    params:
        options = lambda wildcards: get_stat_options(wildcards, rulename="popoolation2_snp_frequency_diff"),
        samples = lambda wildcards: pools if wildcards.sex == "common" else pools[pools.sex.isin([wildcards.sex])]
    threads: lambda wildcards: resources("popoolation2_snp_frequency_diff", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{statname}/{sex}.{region}.{target}.sync_rc.txt{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/snp_frequency_diff"



rule popoolation2_gather_parallel_results:
    """Gather results from parallel analyses"""
    output: analysis = "{results}/{group}/analysis/{analysis}/{statname}/{sex}.{region}{tag}.{suffix}"
    input: unpack(popoolation2_gather_parallel_results_input)
    wildcard_constraints:
        suffix = "(fet.txt.gz|fst.txt.gz|sync_rc.txt.gz|sync_pwc.txt.gz)",
        tag = "(|\.w\d+\.s\d+)"
    resources:
        runtime = resources("popoolation2_gather_parallel_results", "runtime")
    params:
        options = get_params("popoolation2_gather_parallel_results", "options")
    log: "logs/{results}/{group}/analysis/{analysis}/{statname}/{sex}.{region}{tag}.{suffix}.log"
    threads: get_params("popoolation2_gather_parallel_results", "threads")
    shell:
        "cat {input.analysis} > {output.analysis}"
