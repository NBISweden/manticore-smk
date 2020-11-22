def all_rawvc_input(wildcards):
    regions = list(config['workflow']['regions'].keys())
    pfx = str(__RESULTS__ / "ind/rawvc/gatkhc/{region}.vcf.gz")
    val = expand(pfx, region=regions)
    tbi = [f"{x}.tbi" for x in val]
    return {'rawvc': val, 'rawvc.tbi': tbi}


def rawvc_gatkhc_targets_input(wildcards):
    """Get targets for raw variant calling for a sample"""
    ref = config['db']['ref']
    fai = f"{ref}.fai"
    targets = os.path.join(
        os.path.dirname(ref), "gatkhc", f"{wildcards.region}.{wildcards.target}.bed")
    bam = map_dedup_sample_target(wildcards)
    bai = [f"{x}.bai" for x in bam]
    faext = wildcards_or(ext["fa"])
    d = re.sub(faext, ".dict", ref)
    return {'ref': ref, 'fai': fai, 'targets': targets,
            'bam': bam, 'bai': bai, 'dict': d}


def all_gatkhc_samples(wildcards):
    d = dict(wildcards)
    fn = str(__INTERIM__ / "{group}/rawvc/gatkhc/{{SM}}.{target}.{region}.g.vcf.gz").format(**d)
    return expand(fn, SM=individuals["SM"].tolist())


def rawvc_gatk_genomics_db_import_input(wildcards):
    ref = config['db']['ref']
    # Want vcfs for all samples, one region
    vcf = all_gatkhc_samples(wildcards)
    tbi = [f"{x}.tbi" for x in vcf]
    targets = os.path.join(
        os.path.dirname(ref), "gatkhc", f"{wildcards.region}.{wildcards.target}.bed")
    return {'vcf': vcf, 'tbi': tbi, 'targets': targets, 'ref': ref}


def rawvc_gatk_genotype_gvcfs_input(wildcards):
    """Retrieve genomics db databases from raw variant calling"""
    ref = config['db']['ref']
    db = f"{wildcards.results}/{wildcards.group}/rawvc/gatkhc/genomicsdb/{wildcards.region}.{wildcards.target}.db"
    faext = wildcards_or(ext["fa"])
    d = re.sub(faext, ".dict", ref)
    targets = os.path.join(
        os.path.dirname(ref), "gatkhc", f"{wildcards.region}.{wildcards.target}.bed")
    return {'ref': ref, 'db': db, 'dict': d, 'targets': targets}


def _rawvc_vcfs_targets_input(wildcards):
    """Generic function to generate vcf targets to be merged"""
    pfx = f"{str(__RESULTS__)}/{wildcards.group}/rawvc/gatkhc/{wildcards.region}.{{target}}.vcf.gz"
    npart = config['workflow']['regions'].get(wildcards.region, {}).get('npart', 1)
    infiles = expand(pfx, target=list(range(npart)))
    tbi = [f"{x}.tbi" for x in infiles]
    return {'vcfs': infiles, 'tbi': tbi}



def rawvc_picard_merge_vcfs_targets_input(wildcards):
    """Input merge targets for picard. Deprecated in favor of bcftools concat"""
    return _rawvc_vcfs_targets_input(wildcards)


def rawvc_bcftools_concat_vcfs_targets_input(wildcards):
    """Input merge targets for bcftools concat"""
    return _rawvc_vcfs_targets_input(wildcards)
