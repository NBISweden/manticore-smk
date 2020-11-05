bwa_ext = [".amb", ".ann", ".bwt", ".pac", ".sa"]

def bwa_mem_rg(wildcards):
    df = reads[(reads["id"] == 1) & (reads["SM"] == wildcards.sample) & (reads["unit"] == wildcards.unit)]
    rg = {
        'LB': 'lib',
        'PL': 'ILLUMINA',
        'SM': wildcards.sample,
        'PU': wildcards.unit
    }
    rglist = df.to_dict('rglist')
    assert len(rglist) == 1, "only one sample and unit should match read configuration"
    rg.update(rglist[0])
    rg['id'] = config['reads']['read_group_fmt'].format(**rglist[0])
    rgstr = r'-R "@RG\tID:{id}LB:{LB}\tPL:{PL}\tSM:{SM}\tPU:{PU}"'.format(**rg)
    return rgstr


def bwa_mem_index(wildcards):
    """Retrieve bma mem index"""
    genome = os.path.splitext(os.path.basename(config['db']['ref']))[0]
    return f"{wildcards.interim}/map/bwa/index/{genome}"


def bwa_mem_index_ext(wildcards):
    """Retrieve bwa mem index with extensions"""
    index = bwa_mem_index(wildcards)
    return expand("{index}{ext}", index=index, ext=bwa_ext)


def bwa_mem_input(wildcards):
    df = reads[(reads["SM"] == wildcards.sample) & (reads["unit"] == wildcards.unit)]
    return sorted(df['reads.trimmed'].to_list())


def bwa_mem_sample(wildcards):
    df = reads[(reads["id"] == 1) & (reads["SM"] == wildcards.sample)]
    fn = str(__INTERIM__ / "map/bwa/{SM}.bam")
    bam = [fn.format(**x) for k, x in df.iterrows()]
    return bam


def all_bwa_mem_samples(wildcards):
    df = reads[reads["id"] == 1]
    fn = str(__INTERIM__/ "map/bwa/{SM}.bam")
    bam = [fn.format(**x) for k, x in df.iterrows()]
    return bam


def picard_merge_sam_input(wildcards):
    df = reads[(reads["SM"] == wildcards.sample) & (reads["id"] == 1)]
    fn = str(__INTERIM__/ "map/bwa/{SM}/{unit}.bam")
    bam = [fn.format(**x) for k, x in df.iterrows()]
    return bam
