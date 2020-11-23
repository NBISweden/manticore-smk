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
    output: sync = "{results_pool}/raw/popoolation2/{sex}.{region}.{target}.sync{gz}"
    input: mpileup = __INTERIM_POOL__ / "raw/popoolation2/{sex}.{region}.{target}.mpileup.gz"
    resources:
        mem_mb = lambda wildcards, attempt: resources("popoolation2_mpileup2sync_jar", "mem_mb", attempt),
        runtime = lambda wildcards, attempt: resources("popoolation2_mpileup2sync_jar", "runtime", attempt)
    params:
        java_options = get_params("popoolation2_mpileup2sync_jar", "java_options")
    threads: get_params("popoolation2_mpileup2sync_jar", "threads")
    log: "logs/{results_pool}/raw/popoolation2/{sex}.{region}.{target}.sync{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/mpileup2sync_jar"


rule popoolation2_indel_filtering_identify_indel_regions:
    """Identify indel regions in mpileup file"""
    output: gtf = "{results_pool}/raw/popoolation.indels/{sex}.{region}.{target}.mpileup{gz}.indels.gtf"
    input: mpileup = __INTERIM_POOL__ / "raw/popoolation2/{sex}.{region}.{target}.mpileup{gz}"
    resources:
        mem_mb = lambda wildcards, attempt: resources("popoolation2_indel_filtering_identify_indel_regions", "mem_mb", attempt),
        runtime = lambda wildcards, attempt: resources("popoolation2_indel_filtering_identify_indel_regions", "runtime", attempt)
    params:
        java_options = get_params("popoolation2_indel_filtering_identify_indel_regions", "java_options")
    threads: get_params("popoolation2_indel_filtering_identify_indel_regions", "threads")
    log: "logs/{results_pool}/raw/popoolation.indels/{sex}.{region}.{target}.mpileup{gz}.indels.gtf.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/indel_filtering_identify_indel_regions"


rule popoolation2_gather_parallel_results:
    """Gather results from parallel analyses"""
    output: analysis = "{interim_pool}/popoolation2/{region}/{prefix}{repeatmask}.{filters}{analysis}.{suffix}"
    input: unpack(popoolation2_gather_parallel_results_input)
    wildcard_constraints:
        analysis = "(w\d+.s\d+.fst|sync_rc|sync_pwc|w\d+.s\d+.fet|mpileup.gz.indels)",
        prefix = "(all|male/all|female/all)",
        filters = "([\.a-z]+\.|)",
        suffix = "(gz|gtf)"
    resources:
        runtime = resources("popoolation2_gather_parallel_results", "runtime")
    params:
        options = get_params("popoolation2_gather_parallel_results", "options")
    log: "logs/{interim_pool}/popoolation2/{region}/{prefix}{repeatmask}.{filters}{analysis}.{suffix}.log"
    threads: get_params("popoolation2_gather_parallel_results", "threads")
    shell:
        "cat {input.analysis} > {output.analysis}"
