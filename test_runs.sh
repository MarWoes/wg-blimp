#!/bin/bash

set -e
set -o pipefail

NUM_CORES_AVAILABLE=$(nproc --all)

CREATED_CONFIG_TEST_DIR=$(mktemp -d)

wg-blimp create-config \
  --cores-per-job=$NUM_CORES_AVAILABLE \
  $(realpath test/fastq/) \
  $(realpath /scripts/test/lambda-phage.fa) \
  simulated1,simulated2 \
  simulated3,simulated4 \
  $CREATED_CONFIG_TEST_DIR  \
  $CREATED_CONFIG_TEST_DIR/config.yaml

echo "[INFO] Content of created config:"
cat $CREATED_CONFIG_TEST_DIR/config.yaml

echo simulated1 >  $CREATED_CONFIG_TEST_DIR/g1.txt
echo simulated2 >> $CREATED_CONFIG_TEST_DIR/g1.txt
echo simulated3 >  $CREATED_CONFIG_TEST_DIR/g2.txt
echo simulated4 >> $CREATED_CONFIG_TEST_DIR/g2.txt

wg-blimp create-config \
  --cores-per-job=$NUM_CORES_AVAILABLE \
  --use-sample-files \
  $(realpath test/fastq/) \
  $(realpath /scripts/test/lambda-phage.fa) \
  $CREATED_CONFIG_TEST_DIR/g1.txt \
  $CREATED_CONFIG_TEST_DIR/g2.txt \
  $CREATED_CONFIG_TEST_DIR  \
  $CREATED_CONFIG_TEST_DIR/config.yaml

echo "[INFO] Running snakemake from created config"

echo "[INFO] Dry run"
wg-blimp run-snakemake-from-config --cores=$NUM_CORES_AVAILABLE --dry-run $CREATED_CONFIG_TEST_DIR/config.yaml

echo "[INFO] Actual run"

wg-blimp run-snakemake-from-config --cores=$NUM_CORES_AVAILABLE $CREATED_CONFIG_TEST_DIR/config.yaml

echo "[INFO] Deleting output"

wg-blimp delete-all-output --yes $CREATED_CONFIG_TEST_DIR/config.yaml

echo "[INFO] Running standard workflow"

wg-blimp run-snakemake \
  --cores=$NUM_CORES_AVAILABLE \
  --use-sample-files \
  $(realpath test/fastq/) \
  $(realpath /scripts/test/lambda-phage.fa) \
  $CREATED_CONFIG_TEST_DIR/g1.txt \
  $CREATED_CONFIG_TEST_DIR/g2.txt \
  $CREATED_CONFIG_TEST_DIR

echo "[INFO] Deleting files again"

wg-blimp delete-all-output --yes $CREATED_CONFIG_TEST_DIR/config.yaml

echo "[INFO] Removing folder"
rm -r $CREATED_CONFIG_TEST_DIR
