#!/usr/bin/env bash
# written by oleg osipenko
# 10012020

# ----- activate ASFFAST -----
source /usr/local/miniconda3/etc/profile.d/conda.sh  # This must be invoked to use conda activate in Singularity
conda update -y -n base -c defaults conda
conda create -y -n ASFFAST python=3 numpy flask pandas matplotlib biopython openpyxl
source /usr/local/miniconda3/etc/profile.d/conda.sh  # This must be invoked to use conda activate in Singularity
source activate ASFFAST