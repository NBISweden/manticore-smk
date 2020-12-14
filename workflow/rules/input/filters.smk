def get_previous_filtering_step(wildcards, application, input=None, gather=False):
    """Retrieve previous step to base analysis on.

    Either look back one step in a filtering chain and take the
    previous result one count down, or use results based on a manually
    specified input step.

    :param snakemake.io.Wildcards wildcards: snakemake wildcards object
    :param string application: raw caller program application
    :param int index: index of target filter. None for statistics and plots.
    :param string input: manually specified input stage
    """
    if input is not None:
        # Currently unsure how this would be specified
        return input
    analysis = f"analysis/{wildcards.analysis}"
    raw = "rawvc" if wildcards.group == "ind" else "raw"
    raw_application = f"{raw}/{application}"
    count = None
    if "filternum" in wildcards.keys():
        count = int(wildcards.filternum.lstrip("0"))
    target = f".{wildcards.target}" if "target" in wildcards.keys() else ""
    if gather:
        target = ".{target}"# if target == "" else f".{target}"
    sample = f"{wildcards.sample}." if "sample" in wildcards.keys() else ""
    sex = f"{wildcards.sex}." if "sex" in wildcards.keys() else ""
    if count == 1:
        fmt = f"{{root}}/{raw_application}/{sex}{sample}{wildcards.region}{target}{{ext}}"
        return fmt
    nfilters = len(config[analysis]["filters"])
    ## Look at last step if count is none (for stats and plots)
    if count is None:
        count = nfilters + 1
    ## Look at previous filtering step
    num = str(count-1).zfill(2)
    previous_step = list(config[analysis]["filters"].keys())[count - 2]
    previous_filter = config[analysis]["filters"][previous_step]["filter"]
    fmt = f"{{root}}/{analysis}/{num}_{previous_step}_{previous_filter}/{sex}{sample}{wildcards.region}{target}{{ext}}"
    return fmt


def _get_vcf_tbi_input(wildcards, gather=False):
    fmt = get_previous_filtering_step(wildcards, "gatkhc", gather=gather)
    npart = config['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    root = f"{wildcards.results}/{wildcards.group}"
    val = {
        'vcf': expand(fmt, root=root, ext=".vcf.gz", target=list(range(npart))),
        'tbi': expand(fmt, root=root, ext=".vcf.gz.tbi", target=list(range(npart)))
    }
    return val


def filter_bcftools_concat_vcfs_input(wildcards):
    """Get input for filter_bcftools_concat_vcfs"""
    return _get_vcf_tbi_input(wildcards, gather=True)


def filter_gatk_jexl_filter_variants_input(wildcards):
    """Get input for filter_gatk_jexl_filter_variants"""
    return _get_vcf_tbi_input(wildcards)


def filter_gatk_select_variants_input(wildcards):
    """Get input for filter_gatk_select_variants"""
    return _get_vcf_tbi_input(wildcards)


def manticore_plotvcf_input(wildcards):
    """Get input for manticore_plotvcf"""
    return _get_vcf_tbi_input(wildcards)


def get_popoolation_filter_input(wildcards, **kwargs):
    """Get input for popoolation based filtering rules"""
    fmt = get_previous_filtering_step(wildcards, "popoolation")
    npart = config['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    root = f"{wildcards.results}/{wildcards.group}"
    val = {
        'pileup': expand(fmt, root=root, ext=".pileup.gz", target=list(range(npart))),
    }
    return val

def get_popoolation2_filter_input(wildcards, **kwargs):
    """Retrieve partitioned input for popoolation2 filters"""
    fmt = get_previous_filtering_step(wildcards, "popoolation2")
    npart = config['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    root = f"{wildcards.results}/{wildcards.group}"
    val = {
        'sync': expand(fmt, root=root, ext=".sync.gz", target=list(range(npart))),
    }
    return val
