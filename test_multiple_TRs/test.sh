#!/bin/bash

mTR=../mTR

# Synthetic reads with two types of tandem repeats.
# The lengths of two tandem repeats are i and j.
# File named i_j_unit.txt contain the unit strings of length i and j.
for i in 3 5 10 20 50; do
for j in 3 5 10 20 50; do
if [ $i -lt $j ]; then
echo $i"_"$j"-----------------------------------"
$mTR "data/"$i"_"$j".fasta"
fi
done
done

# Reads with four types of tandem repeats
echo "2_5_10_20 -----------------------------------"
$mTR data/2_5_10_20_set.fasta

# Reads with seven types of tandem repeats
echo "2_5_10_20_50_100_200 -----------------------------------"
$mTR data/2_5_10_20_50_100_200_set.fasta

# Nanopore reads with long tandem repeats from C. elegans strain PD1074
echo "-----------------------------------"
$mTR data/worm_chrI.fasta
echo "-----------------------------------"
$mTR data/worm_chrII_1.fasta
echo "-----------------------------------"
$mTR data/worm_chrII_2.fasta

exit 0
