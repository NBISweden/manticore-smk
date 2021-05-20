import os
import re
import sys
import numpy as np
import pandas as pd
import contextlib
from collections import OrderedDict
from snakemake.utils import logger, validate
from snakemake.io import _load_configfile

WORKFLOW_DIR = workflow.current_basedir
SCHEMA_DIR = os.path.realpath(os.path.join(WORKFLOW_DIR, os.pardir, os.pardir, "schemas"))

def wildcards_or(items, empty=False):
    items = list(set(items))
    if empty:
        items = [""] + items
    return f'({"|".join(items)})'


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


class Schema(OrderedDict):
    def __init__(self, data):
        if isinstance(data, str):
            data = _load_configfile(data, filetype="Schema")
        elif not isinstance(data, OrderedDict):
            raise
        for k, v in data.items():
            if isinstance(v, OrderedDict):
                data[k] = Schema(v)
        super().__init__(data)

    def __getattr__(self, attr):
        if attr in self:
            return self.__getitem__(attr)
        try:
            return super().__getattribute__(attr)
        except AttributeError as e:
            raise e


class FilterSchema(Schema):
    def __init__(self, data):
        super().__init__(data)

    @property
    def names(self):
        values = []
        for k in self.definitions.analysisfilters["items"].oneOf:
            fkey = re.sub(".+/", "", k["$ref"])
            values.extend(self.definitions[fkey].properties.keys())
        return values


sample_schema = Schema(os.path.join(SCHEMA_DIR, "samples.schema.yaml"))
analysis_schema = Schema(os.path.join(SCHEMA_DIR, "analysisset.schema.yaml"))
statistic_schema = Schema(os.path.join(SCHEMA_DIR, "statistic.schema.yaml"))
filter_schema = FilterSchema(os.path.join(SCHEMA_DIR, "filter.schema.yaml"))
definitions = Schema(os.path.join(SCHEMA_DIR, "definitions.schema.yaml"))


class PloidyException(Exception):
    pass

class SampleData:
    _index = ["SM"]
    _schemafile = os.path.join(SCHEMA_DIR, "samples.schema.yaml")

    def __init__(self, *args, ignore=None):
        self._data = pd.DataFrame()
        if len(args) == 1:
            args = args[0]
            if isinstance(args, SampleData):
                self._index = args._index
                self._schemafile = args._schemafile
                self._data = args.data
            elif isinstance(args, str):
                self._read_tsv(args)
        elif len(args) > 1:
            assert all(isinstance(x, SampleData) for x in args), \
                logger.error("all instances must be SampleData")
            self._index = args[0]._index
            self._schemafile = args[0]._schemafile
            self._data = pd.concat(x.data for x in args)
        else:
            raise TypeError
        validate(self.data, schema=self.schemafile)
        if ignore is not None:
            self._data = self.subset(samples=ignore, invert=True).data


    def _read_tsv(self, infile):
        self._data = pd.read_csv(infile, sep="\t").set_index(self._index, drop=False)
        self._data = self._data.replace({np.nan: None})
        self._data.index.names = self._index

    @property
    def schemafile(self):
        return self._schemafile

    @property
    def data(self):
        return self._data

    @property
    def samples(self):
        return self.data.SM

    @property
    def unique_samples(self):
        return self.data.SM[list(set(self.data.SM.tolist()))]

    @property
    def samplesize(self):
        return self.data["size"]

    @property
    def populations(self):
        return self.data.population

    @property
    def sex(self):
        return self.data.sex

    def subset(self, invert=False, **kw):
        keep = self.data.SM.isin(self.data.SM)
        for k, v in kw.items():
            if isinstance(v, set):
                v = list(v)
            if not isinstance(v, list):
                v = [v]
            if k == "samples":
                keep = keep & self.data.SM.isin(v)
            else:
                try:
                    keep = keep & self.data[k].isin(v)
                except KeyError as e:
                    print(e)
                    raise
        cls = type(self)
        new = cls(self)
        if invert:
            new._data = new._data[~keep]
        else:
            new._data = new._data[keep]
        return new


class ReadData(SampleData):
    _index = ["SM", "unit", "id"]
    _schemafile = os.path.join(SCHEMA_DIR, "reads.schema.yaml")

    def __init__(self, *args, ignore=None, trim=False):
        super().__init__(*args, ignore=ignore)
        self._data["reads.trimmed"] = self._data["reads"]

    def trim(self, regex_from, regex_to):
        if len(self._data) > 0:
            self._data["reads.trimmed"] = self._data["reads"].str.replace(str(regex_from), str(regex_to / "map/trim"))

    @property
    def populations(self):
        raise NotImplementedError

    @property
    def read_pairs(self):
        """List reads annotated as paired-end (pe) or single-end (se)"""
        df = self.data.groupby(level=["SM", "unit"]).size().to_frame("pe")
        df = self.data[self.data["id"] == 1].droplevel("id").join(df["pe"])
        df["pe"] = df["pe"].map({1: "se", 2: "pe"})
        return df


class Filter(dict):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        assert len(self.keys()) == 1,  f"more than one key defined: {v}"
        self._filter = list(self.keys())[0]
        if self[self._filter] is None:
            self[self._filter] = OrderedDict({})

    @property
    def filter(self):
        return self._filter

    def get(self, attr):
        """Get an attribute from the default key"""
        return self[self.filter].get(attr, None)



class Statistic(dict):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)


class Plot(dict):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)


class Analysis(dict):
    _section = "analysis"

    def __init__(self, name, data, regions, *args, **kw):
        super().__init__(data, *args, **kw)
        self._name = name
        self["filters"] = [Filter(x) for x in self["filters"]]
        self._regions = regions
        self._sex = sample_schema['properties']['sex']['enum']


    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def subset(self, by):
        pass

    @property
    def filters(self):
        return self["filters"]

    @property
    def name(self):
        return self._name.lstrip(self._section + "/")

    @property
    def longname(self):
        return self._name

    def get_filter_options(self):
        pass


class ConfigRule(dict):
    def __init__(self, name, attempt=None, *args, **kw):
        self._name = name
        self._attempt = attempt
        super().__init__(*args, **kw)

    @property
    def attempt(self):
        if self._attempt is None:
            return 1
        return self._attempt

    @property
    def name(self):
        return self._name

    @property
    def threads(self):
        return self.attempt * self["threads"]

    def resources(self, resource):
        assert isinstance(self[resource], int), f"{self}: resource '{resource}' is not an int"
        return self.attempt * self[resource]

    def params(self, attr):
        return self[attr]

    def __str__(self):
        return "ConfigRule: " + self.name


class Config(dict):
    _analysissection = "analysis"

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        for k in self.keys():
            if k.startswith(f"{self._analysissection}/"):
                self[k] = Analysis(k, self[k], regions=self.regions)

    def __getattr__(self, attr):
        if attr in self:
            return self.__getitem__(attr)
        try:
            return super().__getattribute__(attr)
        except AttributeError as e:
            raise e

    def rule(self, rulename, attempt=None):
        """Retrieve rule configuration"""
        ruleobj = ConfigRule(rulename, attempt, self['resources.default'])
        if rulename in self['resources']:
            ruleobj.update(**self['resources'][rulename])
        return ruleobj

    # Region should be a class
    @property
    def regions(self):
        return list(self["workflow"]["regions"].keys())

    def ploidy(self, sample, region, sex=None):
        rconf = self["workflow"]["regions"][region]
        ploidy = rconf["ploidy"].get("common", 2)
        try:
            sex = self["__allsamples__"].data.sex.at[sample]
        except KeyError as e:
            logger.error("No such sample: %s", e)
        try:
            ploidy = rconf["ploidy"][sex]
        except KeyError as e:
            logger.error("No such sex: %s", e)
        finally:
            logger.info(f"falling back on common ploidy {ploidy}")
        return ploidy

    # This should be obtained from FilterSchema?
    def variant_filters(self, rule, vartype):
        filters = self["resources"][rule]["filters"][vartype]
        d = {f"'GATKStandard({v.split()[0]})'": v for v in filters}
        return d

    def get_analysis(self, name):
        return self[f"{self._analysissection}/{name}"]

    @property
    def analyses(self):
        """Retrieve analyses"""
        return [self[k] for k in self.keys() if k.startswith(f"{self._analysissection}/")]

    @property
    def analysisnames(self):
        """Retrieve a list of all analysis names defined in config file"""
        return [x.name for x in self.analyses]

    # Reimplement input/map.smk:bwa_mem_rg
    def read_group(wildcards):
        pass

    @property
    def genome(self):
        """Return short genome name"""
        return os.path.splitext(os.path.basename(self.db["ref"]))[0]



def get_filter_options(wildcards, key="options"):
    """Get filter options."""
    analysis = f"analysis/{wildcards.analysis}"
    index = int(wildcards.filternum.lstrip("0")) - 1
    currentfilter = config[analysis]["filters"][index][wildcards.filtername]
    val = currentfilter.get(key, "")
    return val


def get_filter_input(wildcards):
    """Get filter input."""
    analysis = f"analysis/{wildcards.analysis}"
    index = int(wildcards.filternum.lstrip("0")) - 1
    currentfilter = config[analysis]["filters"][index][wildcards.filtername]
    return currentfilter.get("input", {})


def get_stat_options(wildcards, rulename=None, key="options"):
    """Get stat tool options"""
    options = ""
    if rulename is not None:
        options = get_params(rulename, "options")
    analysis = f"analysis/{wildcards.analysis}"
    index = int(wildcards.statnum.lstrip("0")) - 1
    statistics = config[analysis]["statistics"][index][wildcards.statname]
    val = statistics.get(key, options)
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


# def analysis_subset_regions(key):
#     """Subset regions for a given analysis"""
#     allregions = list(config["workflow"]["regions"].keys())
#     regions = config[key].get("regions", allregions)
#     try:
#         assert set(regions) <= set(allregions)
#     except AssertionError:
#         logger.error(f"configuration section '{key}': some regions undefined: '{regions}'")
#         raise
#     return regions


# def analysis_subset_sex(key, df):
#     """Subset sex for a given analysis"""
#     allsex = df["sex"].tolist() + ["common"]
#     sex = [config[key].get("sex", "common")]
#     try:
#         assert set(sex) <= set(allsex)
#     except AssertionError:
#         logger.error(f"configuration section '{key}': some sexes undefined: '{sex}'")
#         raise
#     return sex


# def analysis_subset_samples(key, df):
#     """Subset samples for a given analysis based on samples and sex keys"""
#     allsamples = df.samples
#     samplelist = config[key].get("samples", allsamples)
#     sex = config[key].get("sex", None)
#     try:
#         assert set(samplelist) <= set(allsamples)
#     except AssertionError:
#         logger.error(f"configuration section '{key}': some samples undefined: '{samplelist}'")
#         raise
#     new = df.subset(samplelist, sex)
#     return new


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
