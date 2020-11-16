import contextlib

class PloidyException(Exception):
    pass


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
