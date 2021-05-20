def all_rawvc_input(wildcards):
    val = {}
    if len(individuals.data) == 0:
        return val
    regions = list(cfg['workflow']['regions'].keys())
    pfx = str(__RESULTS__ / "ind/rawvc/gatkhc/{region}.vcf.gz")
    val = expand(pfx, region=regions)
    pfx = str(__RESULTS__ / "ind/rawvc/gatkhc/{population}.{region}.vcf.gz")
    val = expand(pfx, region=regions, population=individuals.populations)
    tbi = [f"{x}.tbi" for x in val]
    return {'rawvc': val, 'rawvc.tbi': tbi}


def rawvc_gatkhc_targets_input(wildcards):
    """Get targets for raw variant calling for a sample"""
    ref = cfg['db']['ref']
    fai = f"{ref}.fai"
    targets = os.path.join(
        os.path.dirname(ref), "gatkhc", f"{wildcards.region}.{wildcards.target}.bed")
    bam = map_dedup_sample_target(wildcards)
    bai = [f"{x}.bai" for x in bam]
    faext = wildcards_or(ext["fa"])
    d = re.sub(faext, ".dict", ref)
    return {'ref': ref, 'fai': fai, 'targets': targets,
            'bam': bam, 'bai': bai, 'dict': d}


def all_gatkhc_samples(wildcards, population=None):
    d = dict(wildcards)
    fn = str(__INTERIM__ / "{group}/rawvc/gatkhc/{{SM}}.{target}.{region}.g.vcf.gz").format(**d)
    if population is not None:
        return expand(fn, SM=individuals.subset(population=population).unique_samples)
    return expand(fn, SM=individuals.unique_samples)


def rawvc_gatk_genomics_db_import_input(wildcards):
    ref = cfg['db']['ref']
    population = None if wildcards.population == '' else set([wildcards.population])
    # Want vcfs for all samples, one region, and possibly one population
    vcf = all_gatkhc_samples(wildcards, population=population)
    tbi = [f"{x}.tbi" for x in vcf]
    targets = os.path.join(
        os.path.dirname(ref), "gatkhc", f"{wildcards.region}.{wildcards.target}.bed")
    return {'vcf': vcf, 'tbi': tbi, 'targets': targets, 'ref': ref}


def rawvc_gatk_genotype_gvcfs_input(wildcards):
    """Retrieve genomics db databases from raw variant calling"""
    ref = cfg['db']['ref']
    db = "{results}/{group}/rawvc/gatkhc/genomicsdb/{population}{dot}{region}.{target}.db".format(**dict(wildcards))
    faext = wildcards_or(ext["fa"])
    d = re.sub(faext, ".dict", ref)
    targets = os.path.join(
        os.path.dirname(ref), "gatkhc", f"{wildcards.region}.{wildcards.target}.bed")
    return {'ref': ref, 'db': db, 'dict': d, 'targets': targets}


def _rawvc_vcfs_targets_input(wildcards):
    """Generic function to generate vcf targets to be merged"""
    pfx = f"{wildcards.results}/{wildcards.group}/rawvc/gatkhc/{wildcards.population}{wildcards.dot}{wildcards.region}.{{target}}.vcf.gz"
    npart = cfg['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    infiles = expand(pfx, target=list(range(npart)))
    tbi = [f"{x}.tbi" for x in infiles]
    return {'vcf': infiles, 'tbi': tbi}



def rawvc_picard_merge_vcfs_targets_input(wildcards):
    """Input merge targets for picard. Deprecated in favor of bcftools concat"""
    return _rawvc_vcfs_targets_input(wildcards)


def rawvc_bcftools_concat_vcfs_targets_input(wildcards):
    """Input merge targets for bcftools concat"""
    return _rawvc_vcfs_targets_input(wildcards)
