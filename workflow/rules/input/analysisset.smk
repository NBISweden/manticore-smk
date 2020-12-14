## This is hacky!
engine2output_mapping = {
    'popoolation': "{results}/{group}{analysis}/{{name}}/{{sample}}.{{region}}.w{{{{window_size}}}}.s{{{{step_size}}}}.{{statistic}}.txt.gz",
    'popoolation2': "{results}/{group}{analysis}/{{name}}/{{sex}}.{{region}}.{{statistic}}.txt.gz",
    'popoolation2.window': "{results}/{group}{analysis}/{{name}}/{{sex}}.{{region}}.w{{{{window_size}}}}.s{{{{step_size}}}}.{{statistic}}.txt.gz",
}


plot2output_mapping = {
    'manticore-plotvcf.py': "{results}/{group}{analysis}/plots/{{label}}/plot.{{plottype}}.{{region}}.{{ext}}"
}





def all_analysisset_input(wildcards):

    def _update_wildcards(analysis, group):
        d = {}
        if isinstance(wildcards, snakemake.io.Wildcards):
            d = dict(wildcards)
        if "analysis" not in wildcards:
            d.update(**{'analysis': analysis})
        if "results" not in wildcards:
            d.update(**{'results': __RESULTS__})
        if "group" not in wildcards:
            if group is None:
                group = ""
            else:
                group = f"{group}/"
            d.update(**{'group': group})
        wc = snakemake.io.Wildcards(fromdict=d)
        return wc

    allsamples = pd.concat([pools, individuals])
    val = {}

    for k in config.keys():
        if not k.startswith("analysis/"):
            continue

        group = config[k].get("group", None)
        if group is None:
            df = allsamples
        elif group == "ind":
            df = individuals
        elif group == "pool":
            df = pools

        wc = _update_wildcards(k, group)
        regions = analysis_subset_regions(k)
        df = analysis_subset_samples(k, df)
        val[k] = []
        statistics = config[k].get("statistics", [])
        plots = config[k].get("plots", [])
        if len(statistics) == 0 and len(plots) == 0:
            logger.info(f"No statistics or plots section for {k}: skipping")
            continue
        for key, stat in statistics.items():
            pfx = engine2output_mapping[stat['engine']]
            pfx = pfx.format(**wc)
            res = expand(pfx, sample=df["SM"], name=key, region=regions, statistic=stat['statistic'])
            window_size, step_size = analysis_get_window_config(key, stat)
            if window_size is not None:
                res = expand(res, zip, window_size=window_size, step_size=step_size)
            val[k].extend(res)
        for key, plot in plots.items():
            pfx = plot2output_mapping[plot['script']]
            pfx = pfx.format(**wc)
            ext = plot.get("ext", ['png'])
            res = expand(pfx, sample=df["SM"], region=regions, plottype=plot['type'], ext=ext, label=key)
            val[k].extend(res)
    return val
