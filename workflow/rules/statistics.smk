rule statistic_windowed_vcf:
    """Generic rule to calculate plain statistics on individual samples"""
    output:
        txt = "{results}/{group}/analysis/{analysis}/{statnum}_{statname}_{tool}/{region}.{target}.{statistic}.txt.gz"
    input:
        unpack(statistic_vcf_input)
    wildcard_constraints:
        statname = "windowed"
    params:
        options = lambda wildcards: get_stat_options(wildcards)
    resources:
        runtime = lambda wildcards, attempt: resources("statistic_windowed_vcf", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("statistic_windowed_vcf", "mem_mb", attempt),
    threads: get_params("statistic_windowed_vcf", "threads")
    log: "logs/{results}/{group}/analysis/{analysis}/{statnum}_{statname}_{tool}/{region}.{target}.vcf.gz.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/statistic/windowed_vcf"
