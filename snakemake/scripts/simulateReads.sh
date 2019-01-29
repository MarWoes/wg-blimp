#!/bin/bash

SHERMAN=Sherman_v0.1.7/Sherman
SIMULATION_DIR=simulation
DMR_REFERENCE=Sherman_v0.1.7/dmr-ref/
FLANKING_REFERENCE=Sherman_v0.1.7/flanking-ref/

READ_LENGTH=150
NUMBER_OF_SEQS=50000

rm -rf $SIMULATION_DIR
mkdir -p $SIMULATION_DIR

# simulation of flanking regions where each sample has low methylation
for i in `seq 1 4`; do

  $SHERMAN -cr 5 -l $READ_LENGTH -n $NUMBER_OF_SEQS --genome_folder $FLANKING_REFERENCE -pe

  mv simulated_1.fastq $SIMULATION_DIR/simulated"$i"_1.fastq
  mv simulated_2.fastq $SIMULATION_DIR/simulated"$i"_2.fastq
done

# simulation of dmr regions, simulated{1,2} have low methylation, simulated{3,4} high methylation
for i in `seq 1 2`; do

  $SHERMAN -CG 90 -CH 5 -l $READ_LENGTH -n $NUMBER_OF_SEQS --genome_folder $DMR_REFERENCE -pe

  cat simulated_1.fastq >> $SIMULATION_DIR/simulated"$i"_1.fastq
  cat simulated_2.fastq >> $SIMULATION_DIR/simulated"$i"_2.fastq
done

for i in `seq 3 4`; do

  $SHERMAN -CG 5 -CH 5 -l $READ_LENGTH -n $NUMBER_OF_SEQS --genome_folder $DMR_REFERENCE -pe

  cat simulated_1.fastq >> $SIMULATION_DIR/simulated"$i"_1.fastq
  cat simulated_2.fastq >> $SIMULATION_DIR/simulated"$i"_2.fastq
done

# create a little bit more complex structure

mkdir -p $SIMULATION_DIR/more-data
cp $SIMULATION_DIR/simulated1* $SIMULATION_DIR/more-data

rm -f simulated_1.fastq simulated_2.fastq
