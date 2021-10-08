import os
import re
import sys
import json
import numpy as np
import pandas as pd
import contextlib
from collections import OrderedDict
from snakemake.utils import logger, validate
from snakemake.io import _load_configfile
from typing import Union, List, Dict
from dataclasses import dataclass, field, asdict

WORKFLOW_DIR = str(workflow.current_basedir)
SCHEMA_DIR = os.path.realpath(
    os.path.join(WORKFLOW_DIR, os.pardir, os.pardir, "schemas")
)


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


class PropertyDict(OrderedDict):
    """Simple class that allows for property access"""

    def __init__(self, data=dict()):
        super().__init__(data)
        for k, v in data.items():
            if isinstance(v, dict):
                v = PropertyDict(v)
            elif isinstance(v, list):
                val = []
                for x in v:
                    if isinstance(x, PropertyDict):
                        val.append(x)
                    elif isinstance(x, dict):
                        val.append(PropertyDict(x))
                    else:
                        val.append(x)
                v = val
            else:
                pass
            self[k] = v

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        if key not in dir(dict()):
            try:
                setattr(self, key, value)
            except Exception as e:
                print(e)
                print(key, value)
                raise


class Schema(PropertyDict):
    def __init__(self, schemafile):
        self._schemafile = schemafile
        data = _load_configfile(self.schemafile, filetype="Schema")
        super().__init__(data)

    @property
    def schemafile(self):
        return self._schemafile


class FilterSchema(Schema):
    def __init__(self, schemafile):
        super().__init__(schemafile)

    @property
    def names(self):
        values = []
        for k in self.definitions.analysisfilters["items"].oneOf:
            fkey = re.sub(".+/", "", k["$ref"])
            values.extend(self.definitions[fkey].properties.keys())
        return values


sample_schema = Schema(os.path.join(SCHEMA_DIR, "samples.schema.yaml"))
reads_schema = Schema(os.path.join(SCHEMA_DIR, "reads.schema.yaml"))
analysis_schema = Schema(os.path.join(SCHEMA_DIR, "analysisset.schema.yaml"))
statistic_schema = Schema(os.path.join(SCHEMA_DIR, "statistic.schema.yaml"))
filter_schema = FilterSchema(os.path.join(SCHEMA_DIR, "filter.schema.yaml"))
definitions = Schema(os.path.join(SCHEMA_DIR, "definitions.schema.yaml"))


class PloidyException(Exception):
    pass


class SampleData:
    _index = ["SM"]
    _schemafile = sample_schema.schemafile

    def __init__(self, *args, ignore=None):
        if ignore is None:
            ignore = []
        self._data = pd.DataFrame(columns=self._index)
        if len(args) == 1:
            args = args[0]
            if isinstance(args, SampleData):
                self._index = args._index
                self._schemafile = args._schemafile
                self._data = args.data
            elif isinstance(args, str):
                self._read_tsv(args)
        elif len(args) > 1:
            assert all(isinstance(x, SampleData) for x in args), logger.error(
                "all instances must be SampleData"
            )
            self._index = args[0]._index
            self._schemafile = args[0]._schemafile
            self._data = pd.concat(x.data for x in args)
        else:
            raise TypeError
        validate(self.data, schema=self.schemafile)
        if len(ignore) > 0:
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
            if k == "sex" and v == "common":
                continue
            if k == "group":
                if v == "ind":
                    keep = keep & (self.samplesize == 1)
                elif v == "pool":
                    keep = keep & (self.samplesize != 1)
                else:
                    logger.error("No such group {}".format(v))
                    sys.exit()
                continue
            if isinstance(v, set):
                v = list(v)
            if not isinstance(v, list):
                v = [v]
            if len(v) == 0:
                continue
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
    _schemafile = reads_schema.schemafile

    def __init__(self, *args, ignore=None):
        super().__init__(*args, ignore=ignore)
        self._data["reads.trimmed"] = self._data.get("reads")

    def trim(self, regex_from, regex_to):
        if len(self._data) > 0:
            self._data["reads.trimmed"] = self._data["reads"].str.replace(
                str(regex_from), str(regex_to / "map/trim")
            )

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


##############################
# Config related classes
##############################
class AnalysisItem(PropertyDict):
    _section = None
    _label = None

    def __init__(self, data, analysis, index=0, tool=None):
        if data is None:
            data = PropertyDict(dict(raw=PropertyDict()))
        for k, v in data.items():
            if "tool" not in v.keys():
                v["tool"] = tool
                data[k] = v
        super().__init__(data)
        assert len(self.keys()) == 1, f"only one key allowed: {self.keys()}"
        self._name = list(self.keys())[0]
        self._index = int(index)
        self._analysis = analysis
        self._wildcards = None
        self._scatter = False

    @property
    def wildcards(self):
        if self._wildcards is None:
            return dict()
        d = dict(self._wildcards)
        if d != {} and "target" not in d.keys():
            try:
                npart = self.get_region(d["region"]).npart
            except Exception as e:
                print(e)
                npart = self.npart[0]
            d["target"] = range(npart)
        return d

    @property
    def name(self):
        return self._name

    @property
    def index(self):
        return self._index

    @property
    def num(self):
        return str(self.index).zfill(2)

    @property
    def label(self):
        return self._label

    @property
    def group(self):
        return self._analysis.group

    @property
    def sex(self):
        return self._analysis.sex

    @property
    def tool(self):
        return self.get("tool")

    @property
    def prefix(self):
        if self.index >= 1:
            return os.path.join(
                self._analysis.prefix,
                "_".join([f"{self.label}{self.num}", self.name, self.tool]),
            )
        else:
            # Ugly hack; the general tool gatk is different from the
            # variant caller application name gatkhc
            tool = self.tool
            if tool == "gatk":
                tool = "gatkhc"
            if self.group == "ind":
                return os.path.join(self.results, "rawvc", tool)
            else:
                return os.path.join(self.results, "raw", tool)

    @property
    def results(self):
        return self._analysis.results

    @property
    def populations(self):
        return self._analysis.populations

    @property
    def unique_samples(self):
        return self._analysis.allsamples.unique_samples

    @property
    def regions(self):
        return self._analysis.regions

    @property
    def region_names(self):
        return [r.name for r in self.regions]

    def get_region(self, name):
        return self._analysis.get_region(name)

    @property
    def npart(self):
        return [r.npart for r in self.regions]

    @property
    def scatter(self):
        return self._scatter

    @property
    def previous(self):
        """Return previous step"""
        if self.index >= 1:
            return self._analysis.get(self._section)[self.index - 1]
        return None

    def get(self, attr, default=None):
        """Get an attribute from the default key"""
        return self[self.name].get(attr, default)


class Filter(AnalysisItem):
    _section = "filters"
    _label = "f"

    def __init__(self, data, analysis, index=0, tool=None):
        super().__init__(data, analysis, index=index, tool=tool)

    ## FIXME: hacky solution. Want to return target extension
    ## depending on group and tool. IOW we want a Tool class...
    @property
    def ext(self):
        """File extension for filter results"""
        if self.group == "ind":
            return ".vcf.gz"
        elif self.group == "pool":
            if self.tool == "popoolation2":
                return ".sync.gz"
            elif self.tool == "popoolation":
                return ".pileup.gz"
        else:
            pass

    @property
    def fmt(self):
        target = ""
        if self.scatter:
            target = ".{target}"
        else:
            if "target" in self.wildcards.keys():
                target = f".{self.wildcards['target']}"
        if self.group == "pool":
            self._scatter = True
            if self.tool == "popoolation":
                out = f"{{sample}}.{{region}}{target}{self.ext}"
            elif self.tool == "popoolation2":
                out = f"{{sex}}.{{region}}{target}{self.ext}"
            else:
                raise Exception
        elif self.group == "ind":
            out = f"{{population}}{{dot}}{{region}}{target}{self.ext}"
        else:
            raise Exception
        return os.path.join(self.prefix, out)

    def expand(self, wildcards, previous=True):
        """Expand targets taking wildcards context into consideration"""
        filt = self
        if previous:
            filt = self.previous
            if self.name in ["concat"]:
                filt._scatter = True
        if isinstance(wildcards, snakemake.io.Wildcards):
            filt._wildcards = wildcards
        if filt.wildcards:
            return expand(filt.fmt, **filt.wildcards)
        val = []
        pops = filt.populations + [""]
        dots = ["."] * len(filt.populations) + [""]
        for r in self.regions:
            for pop, dot in zip(pops, dots):
                d = {
                    "population": pop,
                    "dot": dot,
                    "sex": filt.sex,
                    "region": r.name,
                    "sample": filt.unique_samples,
                    "target": list(range(r.npart)),
                }
                val.extend(expand(filt.fmt, **d))
        return val


class Statistic(AnalysisItem):
    _section = "statistics"
    _label = "s"

    def __init__(self, data, analysis, index=0, tool=None):
        super().__init__(data, analysis, index=index, tool=tool)

    @property
    def window_config(self):
        window_size = self.get("window_size")
        step_size = self.get("step_size", window_size)
        try:
            assert len(step_size) == len(window_size)
        except AssertionError:
            logger.error(
                f"config section 'statistics:{self.name}' window size and step size must be of equal lengths"
            )
            raise
        return window_size, step_size

    @property
    def ext(self):
        """File extension for statistics results"""
        return ".txt.gz"

    @property
    def fmt(self):
        fmt = "{region}"
        if self.tool == "popoolation":
            fmt = "{sample}.{region}"
        elif self.tool == "popoolation2":
            fmt = "{sex}.{region}"
        if self.name == "windowed_statistic":
            fmt = f"{fmt}.w{{window_size}}.s{{step_size}}.{{statistic}}"
        else:
            if self.tool == "popoolation2":
                fmt = f"{fmt}.sync_{{statistic}}"
        fmt = f"{fmt}{self.ext}"
        return os.path.join(self.prefix, fmt)

    def expand(self, wildcards, previous=True):
        """Expand targets taking wildcards context into consideration"""
        if isinstance(wildcards, snakemake.io.Wildcards):
            self._wildcards = wildcards
        if self.wildcards:
            return expand(self.fmt, **self.wildcards)
        val = []
        for r in self.regions:
            for statistic in self.get("statistic"):
                d = {
                    "statistic": statistic,
                    "sex": self.sex,
                    "region": r.name,
                    "sample": self.unique_samples,
                }
                if self.name == "windowed_statistic":
                    window_size = self.window_config[0]
                    step_size = self.window_config[1]
                    for w, s in zip(window_size, step_size):
                        d["window_size"] = w
                        d["step_size"] = s
                        val.extend(expand(self.fmt, **d))
                else:
                    val.extend(expand(self.fmt, **d))
        return val


class Plot(AnalysisItem):
    _section = "plots"
    _label = "p"

    def __init__(self, data, analysis, index=0, tool=None):
        super().__init__(data, analysis, index=index, tool=tool)


class Analysis(PropertyDict):
    _section = "analysis"

    def __init__(self, name, default, **kw):
        # Set defaults
        default.update(**kw)
        if default["regions"] == []:
            default["regions"] = default["_regions"]
        else:
            default["regions"] = [
                r for r in default["_regions"] if r.name in default["regions"]
            ]
        super().__init__(default)
        self._results = kw.get("results", "results")
        self._name = name
        # Subset samples already here
        self["allsamples"] = self.allsamples.subset(
            group=self.group, samples=self.samples, sex=self.sex
        )
        # Insert 0th filter and update filters
        self["filters"] = [
            Filter(x, self, index=i + 1, tool=self.tool)
            for i, x in enumerate(self["filters"])
        ]
        self["plots"] = [
            Plot(x, self, index=i + 1, tool=self.tool)
            for i, x in enumerate(self["plots"])
        ]
        self["statistics"] = [
            Statistic(x, self, index=i + 1, tool=self.tool)
            for i, x in enumerate(self["statistics"])
        ]
        self["filters"].insert(0, Filter(None, self, tool=self.tool))
        self["plots"].insert(0, Plot(None, self, tool=self.tool))
        self["statistics"].insert(0, Statistic(None, self, tool=self.tool))

    @property
    def name(self):
        return self._name.lstrip(f"{self._section}/")

    @property
    def longname(self):
        return self._name

    @property
    def populations(self):
        return self.allsamples.data.population.tolist()

    def get_region(self, name):
        for r in self.regions:
            if name == r.name:
                return r
        logger.error(f"No such region '{name}' for analysis {self.name}")
        return None

    @property
    def check(self):
        if (
            len(self.statistics) == 0
            and len(self.plots) == 0
            and len(self.filters) == 0
        ):
            logger.info(
                f"No filters, statistics or plots section for {self.name}: skipping"
            )
            return False
        return True

    @property
    def prefix(self):
        if self.group is None:
            return os.path.join(self._results, self.longname)
        return os.path.join(self._results, self.group, self.longname)

    @property
    def results(self):
        if self.group is None:
            return os.path.join(self._results)
        return os.path.join(self._results, self.group)


@dataclass
class RuleConfig:
    _default = workflow.default_resources.parsed
    name: str
    envmodules: List[str] = field(default_factory=list)
    options: Union[str, list, dict] = ""
    runtime: int = _default.get("runtime", 100)
    threads: int = 1
    attempt: int = 1
    mem_mb: int = _default.get("mem_mb", 1000)
    java_options: str = _default.get("java_options", "")
    java_tmpdir: str = _default.get("java_tmpdir", _default.get("tmpdir", "tmp"))
    window_size: List[int] = field(default_factory=list)
    step_size: List[int] = field(default_factory=list)
    extra: Dict = field(default_factory=dict)

    def xthreads(self, wildcards, attempt):
        return attempt * self.threads

    def xruntime(self, wildcards, attempt):
        return attempt * self.runtime

    def xmem(self, wildcards, attempt):
        return attempt * self.mem_mb

    def params(self, attr):
        return getattr(self, attr)


class Region(PropertyDict):
    def __init__(self, name, *args):
        super().__init__(*args)
        self._name = name

    @property
    def name(self):
        return self._name


class Config(PropertyDict):
    _analysissection = "analysis"

    def __init__(self, conf, samples, *args, **kw):
        super().__init__(conf)
        self["__allsamples__"] = samples
        for k, v in self.workflow.regions.items():
            self.workflow.regions[k] = Region(k, v)

        self._init_analyses()

    def _init_analyses(self):
        for k in self.keys():
            default = {
                "_regions": self.regions,
                "allsamples": self.allsamples,
                "regions": [],
                "sex": "common",
                "group": None,
                "samples": [],
                "statistics": [],
                "filters": [],
                "plots": [],
            }
            if k.startswith(f"{self._analysissection}/"):
                self[k] = Analysis(k, default, **self[k])

    def ruleconf(self, rulename, analysis=None, **kw):
        """Retrieve rule configuration"""
        rckeys = RuleConfig.__dataclass_fields__.keys()
        data = {"name": rulename, **kw}
        if "rules" in self.keys():
            if rulename in self.rules:
                d = self.rules[rulename]
                data.update(**{k: d[k] for k in d.keys() if k in rckeys})
                data["extra"] = {k: d[k] for k in d.keys() if k not in rckeys}
        return RuleConfig(**data)

    @property
    def regions(self):
        return list(self.workflow.regions.values())

    @property
    def allsamples(self):
        return self.get("__allsamples__", None)

    def ploidy(self, region, sample=None, sex=None):
        rconf = self.workflow.regions[region]
        ploidy = rconf["ploidy"].get("common", 2)
        try:
            sex = self.allsamples.data.sex.at[sample]
        except KeyError as e:
            logger.error(f"No such sample {sample}: {e}")
        try:
            ploidy = rconf["ploidy"][sex]
        except KeyError as e:
            logger.error(f"No such sex defined for ploidy in region {rconf.name}: {e}")
        finally:
            logger.info(f"falling back on common ploidy {ploidy}")
        return ploidy

    # This should be obtained from FilterSchema?
    def variant_filters(self, rulename, vartype):
        filters = self.rules[rulename].filters[vartype]
        d = {f"'GATKStandard({v.split()[0]})'": v for v in filters}
        return d

    def get_analysis(self, name):
        return self[f"{self._analysissection}/{name}"]

    @property
    def analyses(self):
        """Retrieve analyses"""
        return [
            self[k] for k in self.keys() if k.startswith(f"{self._analysissection}/")
        ]

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

    ## Return ruleconf for consistency?
    def params(self, wildcards, attr, rulename=None):
        """Get rule parameters and options in analysis context"""
        analysis = self.get_analysis(wildcards.analysis)
        index = int(wildcards.itemnum.lstrip("0"))
        if "filtername" in dict(wildcards).keys():
            analysis_item = analysis.filters[index]
        elif "statname" in dict(wildcards).keys():
            analysis_item = analysis.statistics[index]
        elif "plotname" in dict(wildcards).keys():
            analysis_item = analysis.plots[index]
        default = None
        if rulename is not None:
            default = self.ruleconf(rulename).params(attr)
        return analysis_item.get(attr, default)


## FIXME: move this function to config
def get_filter_input(wildcards):
    """Get filter input."""
    analysis = f"analysis/{wildcards.analysis}"
    index = int(wildcards.itemnum.lstrip("0")) - 1
    currentfilter = config[analysis]["filters"][index][wildcards.filtername]
    return currentfilter.get("input", {})


## FIXME: document why this is needed
## patternProperties is not validated!
def preprocess_config(config):
    """Update analysis section filter defaults"""

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
