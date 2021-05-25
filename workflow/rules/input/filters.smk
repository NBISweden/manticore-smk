def manticore_plotvcf_input(wildcards):
    """Get input for manticore_plotvcf"""
    return _get_vcf_tbi_input(wildcards)

def get_popoolation_filter_input(wildcards, **kwargs):
    """Get input for popoolation based filtering rules"""
    index = -1
    previous = False
    if "filtername" in dict(wildcards).keys():
        previous = True
        index = int(wildcards.itemnum.lstrip("0"))
    filt = cfg.get_analysis(wildcards.analysis).filters[index]
    val = {
        'pileup': filt.expand(wildcards, previous=previous)
    }
    return val

def get_popoolation2_filter_input(wildcards, **kwargs):
    """Retrieve partitioned input for popoolation2 filters"""
    index = -1
    previous = False
    if "filtername" in dict(wildcards).keys():
        previous = True
        index = int(wildcards.itemnum.lstrip("0"))
    filt = cfg.get_analysis(wildcards.analysis).filters[index]
    val = {
        'sync': filt.expand(wildcards, previous=previous)
    }
    return val

def filter_vcf_input(wildcards):
    """Get input for filter vcf generic rules"""
    index = int(wildcards.itemnum.lstrip("0"))
    filt = cfg.get_analysis(wildcards.analysis).filters[index]
    val = {
        'vcf': filt.expand(wildcards),
    }
    val['tbi'] = [f"{x}.tbi" for x in val['vcf']]
    return val
