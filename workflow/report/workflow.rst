manticore-smk workflow results
===========================================================

[manticore-smk](https://github.com/NBISweden/manticore-smk) is a
Snakemake workflow for non model organism variant calling and
population genomics analyses.

Data organization
-----------------

.. code-block:: shell

   project_name/                <- top-level project folder
   |
   ├── config                   <- configuration directory for Snakemake and other things
   │
   ├── data
   │   ├── external             <- data from third party sources
   │   ├── interim              <- Intermediate data that can be safely deleted
   │   ├── metadata             <- metadata describing raw data files
   │   ├── processed            <- Final processed data used for analyses
   │   └── raw                  <- The original immutable data dump to be treated as read-only.
   │
   ├── logs                     <- Collection of log outputs, e.g. from cluster managers
   │
   ├── reports                  <- Generated analyses and articles as html, pdf and more.
   │   └── figures              <- Graphics for use in reports.
   │
   └── results                  <- Final results for sharing with collaborators, typically derived
