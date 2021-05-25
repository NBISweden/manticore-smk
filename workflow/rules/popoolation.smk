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
        npart = lambda wildcards: cfg['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    log: "logs/{prefix}/{pool_vc}/{region}.{partition}.log"
    threads: 1
    wrapper: f"{WRAPPER_PREFIX}/bio/pybedtools/make_bed_targets"


rule popoolation_bedtools_repeatmask:
    """bedtools: apply repeat mask coordinates to interval file"""
    output: "{prefix}.rm.{partition}.bed"
    input: left = "{prefix}.{partition}.bed",
           right = cfg['db']['repeats']
    threads: 1
    log: "logs/{prefix}.rm.{partition}.bed.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/bedtools/subtract"


rule popoolation_samtools_filter_mpileup:
    """Generate filtered samtools mpileup file for a target region for popoolation"""
    output: pileup = "{results_pool}/raw/popoolation/{sample}.{region}.{target}.pileup{gz}"
    input: unpack(popoolation_samtools_filter_mpileup_input)
    resources:
        runtime = lambda wildcards, attempt: cfg.ruleconf("popoolation_samtools_filter_mpileup", attempt).resources("runtime"),
    params:
        cfg.ruleconf("popoolation_samtools_filter_mpileup").params("options")
    threads: cfg.ruleconf("popoolation_samtools_filter_mpileup").threads
    log: "logs/{results_pool}/raw/popoolation/{sample}.{region}.{target}.pileup{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation/samtools_filter_mpileup"


## FIXME: rename to make consistent with generic filters: filter_pileup_mask
rule popoolation_filter_pileup_by_gtf:
    """Filter indels from pileup file by using indel gtf generated by popoolation2 output"""
    output:
        pileup = temp("{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{sample}.{region}.{target}.pileup{gz}")
    input:
        unpack(get_popoolation_filter_input),
        gtf = popoolation_filter_pileup_by_gtf_input
    wildcard_constraints:
        filtername = "mask",
        tool = "popoolation"
    resources:
        runtime = lambda wildcards, attempt: cfg.ruleconf("popoolation_filter_pileup_by_gtf", attempt).resources("runtime")
    params:
        options = lambda wildcards: cfg.params(wildcards, "options", "popoolation_filter_pileup_by_gtf")
    threads: lambda wildcards: cfg.ruleconf("popoolation_filter_pileup_by_gtf").threads
    log: "logs/{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{sample}.{region}.{target}.pileup{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation/filter_pileup_by_gtf"


rule popoolation_subsample_pileup:
    """Subsample pileup file"""
    output:
        pileup = temp("{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{sample}.{region}.{target}.pileup{gz}")
    input:
        unpack(get_popoolation_filter_input)
    wildcard_constraints:
        filtername = "filter",
        tool = "popoolation"
    resources:
        runtime = lambda wildcards, attempt: cfg.ruleconf("popoolation_subsample_pileup", attempt).resources("runtime")
    params:
        options = lambda wildcards: cfg.params(wildcards, "options", "popoolation_subsample_pileup")
    threads: lambda wildcards: cfg.ruleconf("popoolation_subsample_pileup").threads
    log: "logs/{results}/{group}/analysis/{analysis}/f{itemnum}_{filtername}_{tool}/{sample}.{region}.{target}.pileup{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation/subsample_pileup"


rule popoolation_variance_sliding:
    output: txt = temp("{results}/{group}/analysis/{analysis}/s{itemnum}_{statname}_{tool}/{sample}.{region}.w{window_size}.s{step_size}.{measure}.{target}.txt.gz")
    input: unpack(get_popoolation_filter_input)
    wildcard_constraints:
        measure = "(pi|theta|D)",
        tool = "popoolation",
        statname = "windowed_statistic"
    resources:
        runtime = lambda wildcards, attempt: cfg.ruleconf("popoolation_variance_sliding", attempt).resources("runtime")
    params:
        options = lambda wildcards: cfg.params(wildcards, "options", "popoolation_variance_sliding"),
        ploidy = lambda wildcards: cfg.ploidy(wildcards.region, sample=wildcards.sample),
        size = lambda wildcards: pools.samplesize.at[wildcards.sample]
    threads: lambda wildcards: cfg.ruleconf("popoolation_variance_sliding").threads
    log: "logs/{results}/{group}/analysis/{analysis}/s{itemnum}_{statname}_{tool}/{sample}.{region}.w{window_size}.s{step_size}.{measure}.{target}.txt.gz"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation/variance_sliding"


rule popoolation_gather_parallel_results:
    """Gather results from popoolation parallel analyses"""
    output:
        analysis = "{results}/{group}/analysis/{analysis}/s{itemnum}_{statname}_{tool}/{sample}.{region}.w{window_size}.s{step_size}.{measure}.txt.gz"
    input: unpack(popoolation_gather_parallel_results_input)
    wildcard_constraints:
        measure = "(pi|theta|D)",
        filters = "([\.a-z]+\.|)",
        group = "pool"
    resources:
        runtime = lambda wildcards, attempt: cfg.ruleconf("popoolation_gather_parallel_results", attempt).resources("runtime")
    params:
        options = cfg.ruleconf("popoolation_gather_parallel_results").params("options")
    log: "logs/{results}/{group}/analysis/{analysis}/s{itemnum}_{statname}_{tool}/{sample}.{region}.w{window_size}.s{step_size}.{measure}.txt.gz.log"
    threads: cfg.ruleconf("popoolation_gather_parallel_results").threads
    shell:
        "cat {input.analysis} > {output.analysis}"



localrules: popoolation_bedtools_repeatmask
