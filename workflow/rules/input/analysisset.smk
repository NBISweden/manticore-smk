def all_analysisset_input(wildcards):
    val = {}
    for analysis in cfg.analyses:
        if not analysis.check:
            continue
        val[analysis.name] = []
        if len(analysis.filters) > 1 and len(analysis.statistics) == 1 and len(analysis.plots) == 1:
            f = analysis.filters[-1]
            val[analysis.name] = f.expand(wildcards, previous=False)
        if len(analysis.statistics) > 1:
            for stat in analysis.statistics:
                if stat.name == "raw":
                    continue
                val[analysis.name].extend(stat.expand(wildcards, previous=False))
        if len(analysis.plots) > 1:
            for plot in analysis.plots:
                if plot.name == "raw":
                    continue
                val[analysis.name].extend(plot.expand(wildcards, previous=False))
    return val
