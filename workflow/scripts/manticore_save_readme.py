#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import yaml
import pandas as pd
from collections import OrderedDict


README = """
# Manticore readme

# Analysis sets

## Filters

Available filters:

{filters}

"""

wd = snakemake.config["__workflow_basedir__"]
with open(os.path.join(wd, "schemas/analysisset.schema.yaml"), "r") as fh:
    analysisset = yaml.safe_load(fh)
analysisset_df = pd.DataFrame(
    {
        "filter": "bcftools_view",
        "description": [
            o["description"] for k, o in analysisset["definitions"].items()
        ],
    }
)


d = {"filters": analysisset_df.to_markdown()}

inconfig = snakemake.config


def od2dict(d):
    if isinstance(d, dict):
        tmp = dict()
        for k, v in d.items():
            tmp[k] = od2dict(v)
        return tmp
    return d


config = od2dict(inconfig)

with open(snakemake.output[0], "w") as fh:
    fh.write(README.format(**d))
