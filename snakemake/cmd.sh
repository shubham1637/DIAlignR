#!/bin/sh

#NOTE: You will need snakemake and singularity installed.
# For Snakemake see: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
# For Singularity see:https://sylabs.io/guides/3.0/user-guide/installation.html#build-and-install-an-rpm
# Activate Conda
#source ~/anaconda3/etc/profile.d/conda.sh
# Activate conda env with snakemake
#conda activate snakemake

# Snakemake DIAlignR workflow
snakemake --snakefile Snakefile.dialignr --use-singularity -j 4

