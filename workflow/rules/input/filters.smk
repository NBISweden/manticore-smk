def get_filter_input(wildcards, **kwargs):
    """Retrieve partitioned input for filter function"""
    analysis = f"analysis/{wildcards.analysis}"
    val = {'vcf': [], 'tbi': []}
    count = int(wildcards.filternum.lstrip("0"))
    if kwargs.get("input", None) is None and count == 1:
        if wildcards.group == "ind":
            pfx = str(__RESULTS__ / "ind/rawvc/gatkhc/{region}.vcf.gz")
            val['vcf'] = pfx.format(region=wildcards.region)
            val['tbi'] = f"{val['vcf']}.tbi"
        return val
    ## Else look at previous step
    num = str(count).zfill(2)
    step = list(config[analysis]["filters"].keys())[count - 2]
    previous_filter = config[analysis]["filters"][step]["filter"]
    val['vcf'] = f"{wildcards.results}/{wildcards.group}/analysis/{wildcards.analysis}/{num}_{step}/{previous_filter}.{wildcards.region}.vcf.gz"
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
    step = list(config[analysis]["filters"].keys())[num - 1]
    previous_filter = config[analysis]["filters"][step]["filter"]
    if num < 10:
        num = f"0{num}"
    val['vcf'] = f"{wildcards.results}/{wildcards.group}/{analysis}/{num}_{step}/{previous_filter}.{wildcards.region}.vcf.gz"
    val['tbi'] = f"{val['vcf']}.tbi"
    return val
