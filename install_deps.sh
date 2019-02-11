#!/bin/bash

# APT PACKAGES
apt update
apt install -y wget bzip2 git ttf-dejavu libsm6 build-essential

# CONDA INSTALLATION
cd /opt/
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p /opt/miniconda3
rm miniconda.sh

conda config --add channels conda-forge
conda config --add channels bioconda && \
conda env create -f /pipeline/environment.yml

# ALLOW ALL USER IDS TO ACCESS CONDA/R
chmod a+rwx -R /opt/miniconda3
