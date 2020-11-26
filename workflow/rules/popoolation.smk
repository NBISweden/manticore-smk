rule all_popoolation:
    """Run popoolation analyses"""
    input: unpack(all_popoolation_input)


rule all_popoolation_raw:
    """Generate raw popoolation data"""
    input: unpack(all_popoolation_raw_input)


rule popoolation_pybedtools_make_bed_targets:
    output:
        bed = "{prefix}/{pool_vc}/{region}.{partition}.bed"
    input: bed = "{prefix}/{region}.bed"
    params:
        npart = lambda wildcards: config['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    log: "logs/{prefix}/{pool_vc}/{region}.{partition}.log"
    threads: 1
    wrapper: f"{WRAPPER_PREFIX}/bio/pybedtools/make_bed_targets"


rule popoolation_bedtools_repeatmask:
    """bedtools: apply repeat mask coordinates to interval file"""
    output: "{prefix}.rm.{partition}.bed"
    input: left = "{prefix}.{partition}.bed",
           right = config['db']['repeats']
    threads: 1
    log: "logs/{prefix}.rm.{partition}.bed.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/bedtools/subtract"


rule popoolation_samtools_filter_mpileup:
    """Generate filtered samtools mpileup file for a target region for popoolation"""
    output: pileup = "{results_pool}/raw/popoolation/{sample}.{region}.{target}.pileup{gz}"
    input: unpack(popoolation_samtools_filter_mpileup_input)
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation_samtools_filter_mpileup", "runtime", attempt),
    params:
        get_params("popoolation_samtools_filter_mpileup", "options")
    threads: lambda wildcards: resources("popoolation_samtools_filter_mpileup", "threads")
    log: "logs/{results_pool}/raw/popoolation/{sample}.{region}.{target}.pileup{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation/samtools_filter_mpileup"


rule popoolation_filter_pileup_by_gtf:
    """Filter indels from pileup file by using indel gtf generated by popoolation2 output"""
    output:
        pileup = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sample}.{region}.{target}.pileup{gz}")
    input:
        unpack(get_popoolation_filter_input),
        gtf = popoolation_filter_pileup_by_gtf_input
    wildcard_constraints:
        filternum = "[0-9]{2}",
        filtertype = "popoolation_filter_pileup_by_gtf"
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation_filter_pileup_by_gtf", "runtime", attempt)
    params:
        options = lambda wildcards: get_filter_options(wildcards)
    threads: lambda wildcards: resources("popoolation_filter_pileup_by_gtf", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sample}.{region}.{target}.pileup{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation/filter_pileup_by_gtf"


rule popoolation_subsample_pileup:
    """Subsample pileup file"""
    output:
        pileup = temp("{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sample}.{region}.{target}.pileup{gz}")
    input:
        unpack(get_popoolation_filter_input)
    wildcard_constraints:
        filternum = "[0-9]{2}",
        filtertype = "popoolation_subsample_pileup"
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation_subsample_pileup", "runtime", attempt)
    params:
        options = lambda wildcards: get_filter_options(wildcards)
    threads: lambda wildcards: resources("popoolation_subsample_pileup", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{filternum}_{filtername}_{filtertype}/{sample}.{region}.{target}.pileup{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation/subsample_pileup"


rule popoolation_variance_sliding:
    output: txt = temp("{results}/{group}/analysis/{analysis}/{statname}/{sample}.{region}.w{window_size}.s{step_size}.{measure}.{target}.txt.gz")
    input: unpack(get_popoolation_filter_input)
    wildcard_constraints:
        measure = "(pi|theta|D)",
        filters = "([\.a-z]+|)"
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation_variance_sliding", "runtime", attempt)
    params:
        options = lambda wildcards: get_stat_options(wildcards, rulename="popoolation_variance_sliding"),
        samples = pools
    threads: lambda wildcards: resources("popoolation_variance_sliding", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{statname}/{sample}.{region}.w{window_size}.s{step_size}.{measure}.{target}.txt.gz"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation/variance_sliding"


rule popoolation_gather_parallel_results:
    """Gather results from popoolation parallel analyses"""
    output:
        analysis = "{results}/{group}/analysis/{analysis}/{statname}/{sample}.{region}.w{window_size}.s{step_size}.{measure}.txt.gz"
    input: unpack(popoolation_gather_parallel_results_input)
    wildcard_constraints:
        measure = "(pi|theta|D)",
        filters = "([\.a-z]+\.|)",
    resources:
        runtime = resources("popoolation_gather_parallel_results", "runtime")
    params:
        options = get_params("popoolation_gather_parallel_results", "options")
    log: "logs/{results}/{group}/analysis/{analysis}/{statname}/{sample}{region}.w{window_size}.s{step_size}.{measure}.txt.gz.log"
    threads: get_params("popoolation_gather_parallel_results", "threads")
    shell:
        "cat {input.analysis} > {output.analysis}"



localrules: popoolation_bedtools_repeatmask
