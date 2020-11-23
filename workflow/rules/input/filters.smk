def get_filter_input(wildcards, **kwargs):
    """Retrieve partitioned input for filter function"""
    analysis = f"analysis/{wildcards.analysis}"
    val = {'vcf': [], 'tbi': []}
    count = int(wildcards.filternum.lstrip("0"))
    if kwargs.get("input", None) is None and count == 1:
        if wildcards.group == "ind":
            pfx = "{results}/ind/rawvc/gatkhc/{region}.vcf.gz"
            val['vcf'] = pfx.format(region=wildcards.region)
            val['tbi'] = f"{val['vcf']}.tbi"
        return val
    ## Else look at previous step
    num = str(count-1).zfill(2)
    previous_step = list(config[analysis]["filters"].keys())[count - 2]
    previous_filter = config[analysis]["filters"][previous_step]["filter"]
    val['vcf'] = f"{wildcards.results}/{wildcards.group}/analysis/{wildcards.analysis}/{num}_{previous_step}/{previous_filter}.{wildcards.region}.vcf.gz"
    val['tbi'] = f"{val['vcf']}.tbi"
    return val


def get_plot_input(wildcards):
    analysis = f"analysis/{wildcards.analysis}"
    val = {'vcf': [], 'tbi': []}
    if len(config[analysis]["filters"].keys()) == 0:
        if wildcards.group == "ind":
            pfx = str(__RESULTS__ / "ind/rawvc/gatkhc/{region}.vcf.gz")
            val['vcf'] = pfx.format(region=wildcards.region)
            val['tbi'] = f"{val['vcf']}.tbi"
        return val
    num = len(config[analysis]["filters"])
    current_filter = list(config[analysis]["filters"].keys())[num - 1]
    previous_filter = config[analysis]["filters"][current_filter]["filter"]
    num = str(count).zfill(2)
    val['vcf'] = f"{wildcards.results}/{wildcards.group}/{analysis}/{num}_{current_filter}/{previous_filter}.{wildcards.region}.vcf.gz"
    val['tbi'] = f"{val['vcf']}.tbi"
    return val


def get_popoolation_filter_input(wildcards, **kwargs):
    """Retrieve partitioned input for popoolation filters"""
    analysis = f"analysis/{wildcards.analysis}"
    val = {'pileup': []}
    count = int(wildcards.filternum.lstrip("0"))
    if kwargs.get("input", None) is None and count == 1:
        pfx = str(__RESULTS__ / "pool/raw/popoolation/{sample}.{region}.{target}.pileup.gz")
        val['pileup'] = pfx.format(**wildcards)
        return val
    num = str(count-1).zfill(2)
    previous_step = list(config[analysis]["filters"].keys())[count - 2]
    previous_filter = config[analysis]["filters"][previous_step]["filter"]
    val['pileup'] = f"{wildcards.results}/pool/analysis/{wildcards.analysis}/{num}_{previous_step}/{previous_filter}.{wildcards.sample}.{wildcards.region}.{wildcards.target}.pileup.gz"
    return val


def get_popoolation2_filter_input(wildcards, **kwargs):
    """Retrieve partitioned input for popoolation2 filters"""
    pass


def get_popoolation_stat_input(wildcards, **kwargs):
    """Retrieve partitioned input pileup for popoolation stat calculations"""
    analysis = f"analysis/{wildcards.analysis}"
    count = len(config[analysis].get("filters", []))
    val = {'pileup': []}
    if kwargs.get("input", None) is None and count == 0:
        pfx = str(__RESULTS__ / "pool/raw/popoolation/{sample}.{region}.{target}.pileup.gz")
        val['pileup'] = pfx.format(**wildcards)
        return val
    num = str(count).zfill(2)
    previous_step = list(config[analysis]["filters"].keys())[count - 1]
    previous_filter = config[analysis]["filters"][previous_step]["filter"]
    val['pileup'] = f"{wildcards.results}/pool/analysis/{wildcards.analysis}/{num}_{previous_step}/{previous_filter}.{wildcards.sample}.{wildcards.region}.{wildcards.target}.pileup.gz"
    return val
