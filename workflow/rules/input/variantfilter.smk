def all_variantfilter(wildcards):
    val = []
    regions = list(config['workflow']['regions'].keys())

    # Individuals
    if len(individuals) > 0:
        pfx = str(__RESULTS__ / "ind/{callset}/{caller}/filter/{region}.{vartype}.vcf.gz")
        callset = ["rawvc"]
        vartype = ["snp", "indel"]
        caller = list(config['workflow']['variantcallers']['ind'])
        val = val + expand(pfx, callset=callset, caller=caller, region=regions, vartype=vartype)

    # Pools
    if len(pools) > 0:
        caller = list(config['workflow']['variantcallers']['pool'])

    return {'variantfilter': val}
