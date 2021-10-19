# Snakemake workflow: manticore-smk

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build status](https://github.com/NBISweden/manticore-smk/workflows/Tests/badge.svg?branch=main)](https://github.com/NBISweden/manticore-smk/actions?query=workflow%3ATests) ![License](https://img.shields.io/badge/license-MIT-blue.svg)



Snakemake workflow for non model organism variant calling and
population genomics analyses.

If you use this workflow in a paper, don't forget to give credits to
the authors by citing the URL of this repository and, if
available, its DOI (see above; currently N/A).


## Authors

* Per Unneberg (@percyfal)

## Features

Features include but are not limited to:

* Support for several variant callers:
  - GATK4 HaplotypeCaller
  - freebayes
  - bcftools
  - popoolation
* GATK best practice variant calling
* BQSR and VQSR, where knownsites are identified as the intersection
  of call sets from up to three variant callers
* Variant calling in pools
* Genome size estimates
* Preliminary support for sex chromosome identification
* Genome-wide selection scans with scikit-allel
* Schema-verified configuration

Applications and analyses are organized in modules so that
only a subset of analyses can be selected.


## Quickstart


### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a
   template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository)
   the newly created repository to your local system, into the place
   where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files
in the `config/` folder. Adjust `config.yaml` to configure the
workflow execution. In addition, add at least one of the following
files to define your experimental setup:

`reads.tsv`
    File containing read file uris and corresponding metadata, such as
	sample identifier, sequencing unit, and read pair id

`samples.individual.tsv`
	Sample setup for individual sequences, defining samples, populations
	and other metadata

`samples.pool.tsv`
	Sample setup for pool sequences, defining samples, populations, pool
	size and other metadata

`datasources.tsv`
	data-source key-value pairs definining workflow files and sources
	specified as [Uniform Resource Identifiers
	(uri)](https://en.wikipedia.org/wiki/Uniform_Resource_Identifier)

NOTE: the config directory doesn't have to be in the workflow source
directory, in which case snakemake must be invoked with the full path
to the Snakemake file:

	snakemake -s /path/to/manticore-smk/workflow/Snakefile

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

You can also use a [snakemake
profile](https://github.com/snakemake-profiles/) for fine-tuning
executions. For instance, to use the [slurm
profile](https://github.com/Snakemake-Profiles/slurm) run

	cookiecutter https://github.com/Snakemake-Profiles/slurm.git
	snakemake --use-conda --profile slurm --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above. See the [Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/executable.html)
for further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained
interactive HTML report with all results via:

    snakemake --report report.html

The report contains documentation about the workflow.

### Step 6: Contribute back

In case you have also changed or added steps, please consider
contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the
   original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository)
   the fork to your local system, to a different place than where you
   ran your analysis.
3. Copy the modified files from your analysis to the clone of your
   fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not**
   accidentally copy config file contents or sample sheets. Instead,
   manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull
   request](https://help.github.com/en/articles/creating-a-pull-request)
   against the original repository.


## Configuration

For a quick overview of example configuration files, see [config/config.yaml](https://github.com/NBISweden/manticore-smk/blob/main/config/config.yaml) and the test configuration [.test/config/config.yaml](https://github.com/NBISweden/manticore-smk/blob/main/.test/config/config.yaml)

### Schemas

All configuration files are evaluated against [configuration
schemas](https://github.com/NBISweden/manticore-smk/tree/main/workflow/schemas).
The validation ensures configuration keys are populated, which
minimizes configuration overhead for the user. The schemas are
self-documented and define and describe all available configuration
options.

As an example, `workflow/schemas/samples.ind.schema.yaml` defines a
tabular sample input file format for individual samples. There are six
defined properties `SM` (sample), `population`, `species`, `genus`,
`treatment`, and `sex`, of which `SM` and `population` are required.


See the tutorial [understanding
jsonschema](https://json-schema.org/understanding-json-schema/) for an
accessible introduction to schemas.


### Main workflow configuration

The main workflow configuration entries are `db` and `workflow`. `db`
defines various database resources, such as reference sequence (`ref`)
and repeat file (`repeats`). `workflow` configures general workflow
settings, such as what `qc` programs to run, whether to `trim` or not,
and how to parallelize over `regions`. The latter is the most
important subsection. The `regions` configuration keys define region
names, which in turn define a `bed` file listing what regions to
include, `npart` how many partitions the region will be split into
upon parallel processing, and the `ploidy`. This can be utilized to
apply different settings to autosomes versus sex chromosomes:

    workflow:
      regions:
        autosomes:
          bed: data/external/ref/autosomes.bed
          npart: 20
          ploidy: 2
        Y:
          bed: data/external/ref/chrY.bed
          npart: 4
          ploidy: 1

See
[definitions.schema.yaml](https://github.com/NBISweden/manticore-smk/blob/main/workflow/schemas/definitions.schema.yaml)
for definitions.

See also
[config.schema.yaml](https://github.com/NBISweden/manticore-smk/blob/main/workflow/schemas/config.schema.yaml)
for the main configuration sections.

### Rule configuration

Every rule has a corresponding configuration entry under the
configuration key `rules`. Hence, it is possible to modify and amend
program `options` and in some cases other settings. For details you
currently have to inspect the schema files. The following example
shows how to set some properties for `map_bwa_mem`:

    rules:
      map_bwa_mem:
        options: -k 25
        sort: picard
        sort_order: queryname

### Resource configuration

Recent snakemake versions have added options that allow for
fine-tuning of resources via either the command line or
[profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).
Since the workflow consists of many rules, the latter is recommended.
At its simplest, a profile is a directory with a `config.yaml` file.
For more elaborate profile templates for cluster execution, see
[snakemake profiles
documentation](https://github.com/snakemake-profiles/doc).

An example profile configuration could look as follows:

    restart-times: 3
    max-jobs-per-second: 1
    max-status-checks-per-second: 10
    local-cores: 1
    latency-wait: 60
    default-resources:
      - runtime=100
      - mem_mb=6000
      - disk_mb=1000000
    set-threads:
      - map_bwa_mem=20
      - trim_cutadapt_pe=20
      - qc_qualimap_bamqc_pe=10
    set-resources:
      - qc_qualimap_bamqc_pe:mem_mb=12000
      - rawvc_gatkhc_targets:runtime=1200



## Testing

Test cases are in the subfolder `.test`. They are automatically
executed via continuous integration with [Github
Actions](https://github.com/features/actions). To run the tests, cd to
`.test` and issue

	snakemake --use-conda --conda-frontend mamba --show-failed-logs --cores 2 --conda-cleanup-pkgs cache -s ../workflow/Snakefile --wrapper-prefix file://$(pwd)/../workflow/wrappers
