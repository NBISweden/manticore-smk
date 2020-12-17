from collections import OrderedDict
import contextlib

class PloidyException(Exception):
    pass


def preprocess_config(config):
    """Update filter defaults"""
    def _update_section(k, section):
        values = config[k].get(section, [])
        updated = []
        for v in values:
            assert len(v.keys()) == 1, f"more than one key defined in filter: {v}"
            key = list(v.keys())[0]
            if v[key] is None:
                v[key] = OrderedDict({})
            updated.append(v)
        config[k][section] = updated

    for k in config.keys():
        if not k.startswith("analysis/"):
            continue
        _update_section(k, "filters")
        #_update_section(k, "statistics")


def get_app_params(wildcards, rule):
    return config


def read_pairs_dataframe():
    """List reads as paired-end or single-end"""
    df = reads.groupby(level=["SM", "unit"]).size().to_frame("pe")
    df = reads[reads["id"] == 1].droplevel("id").join(df["pe"])
    df["pe"] = df["pe"].map({1: "se", 2: "pe"})
    return df


def resources(rule, resource, attempt=1, wildcards=None, **kwargs):
    """Retrieve resources for rule multiplying the value by attempt"""
    # TODO: Add prior checks to resource
    if config['resources'][rule].get(resource, None):
        val = config['resources'][rule][resource]
    else:
        val = config['resources.default'][resource]
    return attempt * val


def get_ploidy(sample, region, **kwargs):
    """Retrieve ploidy for a given sample and region"""
    d = config["workflow"]["regions"][region]["ploidy"]
    ploidy = d.get("common", 2)
    try:
        if sample in individuals.SM:
            sex = individuals.at[sample, "sex"]
        elif sample in pools.SM:
            sex = pools.at[sample, "sex"]
        if sex in d.keys():
            ploidy = d[sex]
    except PloidyException as e:
        logger.info(e)
        logger.info(f"falling back on common ploidy: {ploidy}")
    return ploidy


def get_params(rule, resource, wildcards=None, **kwargs):
    """Retrieve rule parameters"""
    val = config['resources'][rule].get(resource, None)
    if val is not None:
        return val
    val = config['resources.default'][resource]
    return val


def get_variant_filters(rule, wildcards, **kwargs):
    vartype = wildcards.vartype
    filters = config["resources"][rule]["filters"][vartype]
    d = {f"'GATKStandard({v.split()[0]})'": v for v in filters}
    return d


def get_filternames():
    """Retrieve a list of all filternames defined in config file"""
    filternames = []
    filterdef = load_schema("../schemas/filter.schema.yaml")
    for k in filterdef["definitions"]["analysisfilters"]["items"]["oneOf"]:
        fkey = re.sub(".+/", "", k["$ref"])
        filternames.extend(filterdef["definitions"][fkey]["properties"].keys())
    return list(set(filternames))


def get_analysisnames():
    """Retrieve a list of all analysis names defined in config file"""
    analyses = []
    for key in config.keys():
        if not key.startswith("analysis/"):
            continue
        analyses.append(key.lstrip("analysis/"))
    return analyses


def get_filter_options(wildcards, gatk=False):
    """Get filter options."""
    analysis = f"analysis/{wildcards.analysis}"
    filters = config[analysis]["filters"]
    val = filters[wildcards.filtername].get("options", "")
    if gatk:
        val = {f"'GATKStandard({v.split()[0]})'": v for v in val}
    return val


def get_filter_input(wildcards):
    """Get filter input."""
    analysis = f"analysis/{wildcards.analysis}"
    filters = config[analysis]["filters"]
    return filters[wildcards.filtername].get("input", {})


def get_stat_options(wildcards, rulename=None):
    """Get stat engine options"""
    options = ""
    if rulename is not None:
        options = get_params(rulename, "options")
    analysis = f"analysis/{wildcards.analysis}"
    statistics = config[analysis]["statistics"]
    val = statistics[wildcards.statname].get("options", options)
    return val


def get_plot_options(wildcards, rulename=None):
    """Get stat engine options"""
    options = ""
    if rulename is not None:
        options = get_params(rulename, "options")
    analysis = f"analysis/{wildcards.analysis}"
    plots = config[analysis]["plots"]
    val = plots[wildcards.label].get("options", options)
    return val


def analysis_subset_regions(key):
    """Subset regions for a given analysis"""
    allregions = list(config["workflow"]["regions"].keys())
    regions = config[key].get("regions", allregions)
    try:
        assert set(regions) <= set(allregions)
    except AssertionError:
        logger.error(f"configuration section '{key}': some regions undefined: '{regions}'")
        raise
    return regions


def analysis_subset_sex(key, df):
    """Subset sex for a given analysis"""
    allsex = df["sex"].tolist() + ["common"]
    sex = [config[key].get("sex", "common")]
    try:
        assert set(sex) <= set(allsex)
    except AssertionError:
        logger.error(f"configuration section '{key}': some sexes undefined: '{sex}'")
        raise
    return sex


def analysis_subset_samples(key, df):
    """Subset samples for a given analysis based on samples and sex keys"""
    allsamples = df["SM"].tolist()
    samplelist = config[key].get("samples", allsamples)
    sex = config[key].get("sex", None)
    try:
        assert set(samplelist) <= set(allsamples)
    except AssertionError:
        logger.error(f"configuration section '{key}': some samples undefined: '{samplelist}'")
        raise
    df = df[df["SM"].isin(samplelist)]
    if sex is not None:
        df = df[df["sex"].isin([sex])]
    return df


def analysis_get_window_config(key, stat):
    if "window_size" not in stat.keys():
        return None, None
    window_size = stat.get("window_size", None)
    step_size = stat.get("step_size", window_size)
    try:
        assert len(step_size) == len(window_size)
    except AssertionError:
        logger.error(f"config section '{k}:statistics:{key}' window size and step size must be of equal lengths")
        raise
    return window_size, step_size




# context manager for cd
@contextlib.contextmanager
def cd(path, logger):
    CWD = os.getcwd()
    logger.info("Changing directory from {} to {}".format(CWD, path))

    os.chdir(path)
    try:
        yield
    except:
        logger.warning("Exception caught: ".format(sys.exc_info()[0]))
    finally:
        logger.info("Changing directory back to {}".format(CWD))
        os.chdir(CWD)
