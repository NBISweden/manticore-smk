# Rules
rule popoolation2_samtools_mpileup:
    """Run samtools mpileup for popoolation2"""
    output: pileup = temp("{interim_pool}/popoolation2/{region}/{sex}{ossep}{all}{repeatmask}.{target}.mpileup{gz}")
    input: unpack(popoolation2_samtools_mpileup_input)
    resources:
        runtime = lambda wildcards, attempt: resources("popoolation2_samtools_mpileup", "runtime", attempt)
    params:
        options = get_params("popoolation2_samtools_mpileup", "options")
    threads: lambda wildcards: resources("popoolation2_samtools_mpileup", "threads")
    log: "logs/{interim_pool}/popoolation2/{region}/{sex}{ossep}{all}{repeatmask}.{target}.mpileup{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/samtools_mpileup"


rule popoolation2_mpileup2sync_jar:
    """Convert mpileup output to sync format"""
    output: sync = "{prefix}.sync"
    input: pileup = "{prefix}.mpileup.gz"
    resources:
        mem_mb = lambda wildcards, attempt: resources("popoolation2_mpileup2sync_jar", "mem_mb", attempt),
        runtime = lambda wildcards, attempt: resources("popoolation2_mpileup2sync_jar", "runtime", attempt)
    params:
        java_options = get_params("popoolation2_mpileup2sync_jar", "java_options")
    threads: get_params("popoolation2_mpileup2sync_jar", "threads")
    log: "logs/{prefix}.sync.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/popoolation2/mpileup2sync_jar"
