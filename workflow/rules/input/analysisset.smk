## This is hacky!
plot2output_mapping = {
    'manticore-plotvcf.py': "{results}/{group}{analysis}/plots/{{label}}/plot.{{plottype}}.{{region}}.{{ext}}"
}


def _make_prefix(wc):
    """Constant prefix for all analysis types"""
    if wc.group == "":
        return "{results}/{analysis}".format(**wc)
    return "{results}/{group}/{analysis}".format(**wc)


def _filter2ext(tool, group):
    """File extensions for filters unless specified by filter 'to' keyword"""
    # Future case
    if group == "":
        pass
    if group == "ind":
        return ".vcf.gz"
    if tool == "popoolation":
        return ".pileup.gz"
    return ".sync.gz"

# FIXME: individual outputs should also map to {sex}.{region}
def _filter2fmt(tool, group):
    ext = _filter2ext(tool, group)
    if tool == 'popoolation':
        return f"{{sample}}.{{region}}{ext}"
    if tool == 'popoolation2':
        return f"{{sex}}.{{region}}{ext}"
    return f"{{region}}{ext}"

# Assume .txt.gz extension for all stats
def _stat2fmt(tool, group, window=True):
    fmt = "{region}"
    if tool == "popoolation":
        fmt = "{sample}.{region}"
    elif tool == "popoolation2":
        fmt = "{sex}.{region}"
    if window:
        fmt = f"{fmt}.w{{{{window_size}}}}.s{{{{step_size}}}}.{{statistic}}"
    else:
        if tool == "popoolation2":
            fmt = f"{fmt}.sync_{{statistic}}"
    fmt = f"{fmt}.txt.gz"
    return fmt

def _plot2fmt(tool, group):
    pass


def all_analysisset_input(wildcards):
    val = {}
    for analysis in cfg.analyses:
        if not analysis.check:
            continue
        if len(analysis.filters) > 0:
            f = analysis.filters[-1]
            val[analysis.name] = f.fmt
    return val



## analysissets should probably be mapped to a class!
def all_analysisset_input2(wildcards):
    """Generate inputs from all analysis sets

    Loop analysis sets and generate input names based on:

    1. last filter output (for ind group)
    2. statistics
    3. plots

    """

    def _update_wildcards(analysis, group):
        d = {}
        if isinstance(wildcards, snakemake.io.Wildcards):
            d = dict(wildcards)
        if "analysis" not in wildcards:
            d.update(**{'analysis': analysis})
        if "results" not in wildcards:
            d.update(**{'results': __RESULTS__})
        if "sex" not in wildcards:
            d.update(**{'sex': 'common'})
        if "group" not in wildcards:
            if group is None:
                group = ""
            d.update(**{'group': group})
        wc = snakemake.io.Wildcards(fromdict=d)
        return wc

    # allsamples = pd.concat([pools.data, individuals.data])
    val = {}

    for analysis in cfg.analyses:
        if not analysis.check:
            continue

    for k in config.keys():
        if not k.startswith("analysis/"):
            continue

        group = config[k].get("group", None)
        if group is None:
            df = allsamples.data
        elif group == "ind":
            df = individuals.data
        elif group == "pool":
            df = pools.data

        wc = _update_wildcards(k, group)
        pfx = _make_prefix(wc)
        regions = analysis_subset_regions(k)
        df = analysis_subset_samples(k, df)
        sex = analysis_subset_sex(k, df)
        kw_target = {
            'region': regions,
            'sex': sex,
            'sample': df["SM"]
        }
        #print(kw_target)
        val[k] = []
        filters = config[k].get("filters", [])
        statistics = config[k].get("statistics", [])
        plots = config[k].get("plots", [])
        tool = config[k].get("tool", None)
        if len(statistics) == 0 and len(plots) == 0 and len(filters) == 0:
            logger.info(f"No filters, statistics or plots section for {k}: skipping")
            continue
        # Add last filter output to results if plots and statistics
        # empty for individual analyses only
        if len(statistics) == 0 and len(plots) == 0 and group == "ind":
            num = str(len(filters)).zfill(2)
            lastfilt = filters[-1]
            key = list(lastfilt.keys())[0]
            tool = lastfilt[key].get("tool", tool)
            if tool is None:
                logger.error(f"No tool defined for filter '{key}' in '{k}'; check configuration file")
                sys.exit(1)
            fmt = _filter2fmt(tool, group)
            target_fmt = f"{pfx}/{num}_{key}_{tool}/{fmt}"
            target = expand(target_fmt, **kw_target)
            val[k].extend(target)

        # Assuming there are statistics or plots to be made
        for i in range(len(statistics)):
            stat = statistics[i]
            num = str(i + 1).zfill(2)
            key = list(stat.keys())[0]
            tool = stat.get("tool", tool)
            if tool is None:
                logger.error(f"No tool defined for statistic '{key}' in '{k}'; check configuration file")
                sys.exit(1)
            fmt = _stat2fmt(tool, group, key == "windowed_statistic")
            target_fmt = f"{pfx}/{num}_{key}_{tool}/{fmt}"
            target = expand(target_fmt, statistic=stat[key]['statistic'], **kw_target)
            if key == "windowed_statistic":
                window_size, step_size = analysis_get_window_config(key, stat[key])
                if window_size is not None:
                    target = expand(target, zip, window_size=window_size, step_size=step_size)
            val[k].extend(target)
    return val


        # for key, plot in plots.items():
        #     pfx = plot2output_mapping[plot['script']]
        #     pfx = pfx.format(**wc)
        #     ext = plot.get("ext", ['png'])
        #     res = expand(pfx, sample=df["SM"], region=regions, plottype=plot['type'], ext=ext, label=key)
        #     val[k].extend(res)
