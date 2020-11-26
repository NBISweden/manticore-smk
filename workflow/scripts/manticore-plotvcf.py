#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import shutil
import argparse
import random
import numpy as np
import allel
import zarr
import numcodecs
import matplotlib.pyplot as plt
import seaborn as sns


sns.set_style("white")
sns.set_style("ticks")


populations = ["CHS", "YRI"]
pop_colours = {"CHS": "#FF0000", "YRI": "#008000"}

sample_population = ["CHS", "YRI"]


def do_pca(x, n, ncomp=10):
    vidx = np.random.choice(x.shape[0], n, replace=False)
    vidx.sort()
    y = x.take(vidx, axis=0)
    coords, model = allel.pca(y, n_components=ncomp, scaler="patterson")
    return coords, model


# Taken from http://alimanfoo.github.io/2015/09/28/fast-pca.html
def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in populations:
        flt = sample_population == pop
        ax.plot(
            x[flt],
            y[flt],
            marker="o",
            linestyle=" ",
            color=pop_colours[pop],
            label=pop,
            markersize=6,
            mec="k",
            mew=0.5,
        )
    ax.set_xlabel(
        "PC%s (%.1f%%)" % (pc1 + 1, model.explained_variance_ratio_[pc1] * 100)
    )
    ax.set_ylabel(
        "PC%s (%.1f%%)" % (pc2 + 1, model.explained_variance_ratio_[pc2] * 100)
    )


def fig_pca(coords, model, title, sample_population):
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population)
    # ax = fig.add_subplot(1, 2, 2)
    # plot_pca_coords(coords, model, 2, 3, ax, sample_population)
    ax.legend(bbox_to_anchor=(1, 1), loc="upper left")
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    return fig


if __name__ == "__main__":

    if snakemake.params.options != "":
        options = snakemake.params.options.split(" ")
        sys.argv.extend(options)

    parser = argparse.ArgumentParser(description="manticore-plotvcf option parser")
    parser.add_argument(
        "--subsample-size", metavar="n", type=int, help="subsample size", default=10000
    )
    args = parser.parse_args()

    vcf_path = str(snakemake.input.vcf)
    output = snakemake.output[0]
    outdir = os.path.dirname(output)
    plottype = snakemake.wildcards.type
    dev = snakemake.wildcards.ext

    ## Convert to zarr
    zarr_path = os.path.join(outdir, os.path.basename(vcf_path) + ".zarr")
    allel.vcf_to_zarr(
        vcf_path,
        zarr_path,
        log=sys.stdout,
        fields="*",
        alt_number=8,
        compressor=numcodecs.Blosc(cname="zstd", clevel=1, shuffle=False),
    )

    callset = zarr.open_group(zarr_path, mode="r")
    g = allel.GenotypeChunkedArray(callset["calldata/GT"])
    n = min(len(g), args.subsample_size)
    gn = g.to_n_alt()

    coords, model = do_pca(gn, n, ncomp=2)
    fig = fig_pca(coords, model, "PCA of first four components", sample_population)
    fig.savefig(output)

    shutil.rmtree(zarr_path)
