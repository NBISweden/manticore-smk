#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import allel
import zarr

vcf = snakemake.input.vcf
output = snakemake.output[0]
outdir = os.path.dirname(output)
print(outdir)
zarr_path = os.path.join(outdir, os.path.basename(vcf) + ".zarr")
allel.vcf_to_zarr(vcf, zarr_path, fields="*", log=sys.stdout, overwrite=True)
print("writing ", zarr_path)
