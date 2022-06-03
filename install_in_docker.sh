#!/bin/bash

set -euxo pipefail

apt update
apt install unzip

TEST_DIR=/tmp/testrun/

mkdir -p $TEST_DIR

cd $TEST_DIR

wget https://uni-muenster.sciebo.de/s/7vpqRSEATYcvlnP/download

unzip download

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install --yes --name base wg-blimp r-base==4.1.1

#conda install --yes --name base click h5py pysam r-base==4.1.1 r-data.table r-dt r-ggplot2 r-htmlwidgets r-httpuv r-shiny r-shinydashboard ruamel.yaml 'snakemake-minimal>=5.8' mamba==0.17.0
conda clean --all --yes

#git clone --recursive https://github.com/MarWoes/wg-blimp.git /root/wg-blimp
#pip install /root/wg-blimp

cd $TEST_DIR
wg-blimp create-config --cores-per-job 4 fastq chr22.fasta blood1,blood2 sperm1,sperm2 results config.yaml
wg-blimp run-snakemake-from-config --cores=4 --dry-run config.yaml
wg-blimp run-snakemake-from-config --cores=4 config.yaml

cd /
rm -rf $TEST_DIR
