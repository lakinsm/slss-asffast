#!/usr/bin/env bash
# written by oleg osipenko
# 10012020


# ----- install miniconda -----
cd /opt || exit
wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh
sh Miniconda3-4.6.14-Linux-x86_64.sh -b -p "/usr/local/miniconda3"
