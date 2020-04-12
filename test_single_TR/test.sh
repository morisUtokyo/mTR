#!/bin/bash

# executable modules
mTR=../mTR
mode=-c #-cp
rand_seq=util/rand_seq
count_match=util/count_match
comp_mTR_DP=util/comp_mTR_DP
# temporary directory
tmp=tmp
# error pattern
substitution_rate=1.6 #12.7
insertion_rate=9.0 #3.2
deletion_rate=3.8 #4.7

#substitution_rate=9.7
#insertion_rate=2.9
#deletion_rate=7.5

# number of sequences
num_rand_seqs=1000 #1000 #100
# Length of random sequences before and after the tandem repeat


echo "This program generates a fasta file with "$num_rand_seqs" sequences with tandem repeats according to the following error pattern, predicts the unit string for each sequence, and outputs the number of correct answers in the prediction."

echo "substitution rate = " $substitution_rate "(%)"
echo "insertion rate    = " $insertion_rate "(%)"
echo "deletion rate     = " $deletion_rate "(%)"

# j frequency of units
for j in 10; do
#for j in 10 20 50 100 200; do
# i unit length
for i in 2 5 10 20 50 100 200; do
#for i in 200 200 200; do

rand_seq_len=$(( (i*j) ))
#rand_seq_len=$(( (i*j)/2 ))
# if [ $((i*j)) -lt 1000 ]; then rand_seq_len=$((i*j)) else rand_seq_len=1000 fi


echo "------------------"
echo "Freq units = "$j", Len. units = "$i
echo "Len randam seqs = "$rand_seq_len
$rand_seq $tmp"/rand_seq.fasta" $tmp"/rand_unit.txt" $i $j $substitution_rate $insertion_rate $deletion_rate $rand_seq_len $rand_seq_len $num_rand_seqs
$mTR $mode $tmp"/rand_seq.fasta" > $tmp"/rand_TRs.txt"

echo "Number of correct answers for" $num_rand_seqs "sequences with tandem repeats"
$count_match $tmp"/rand_TRs.txt" $tmp"/rand_unit.txt"
$comp_mTR_DP $tmp"/rand_unit.txt" $tmp"/rand_TRs.txt" $tmp"/comp_mTR_DP_result.txt"
echo "Number of answers with match ratios >= 1"
cat $tmp"/comp_mTR_DP_result.txt"  | awk '$1>=1' | wc -l
echo "Number of answers with match ratios >= 0.99"
cat $tmp"/comp_mTR_DP_result.txt"  | awk '$1>=0.99' | wc -l
echo "Number of answers with match ratios >= 0.98"
cat $tmp"/comp_mTR_DP_result.txt"  | awk '$1>=0.98' | wc -l
echo "Number of answers with match ratios >= 0.96"
cat $tmp"/comp_mTR_DP_result.txt"  | awk '$1>=0.96' | wc -l
echo "Number of answers with match ratios >= 0.94"
cat $tmp"/comp_mTR_DP_result.txt"  | awk '$1>=0.94' | wc -l
done
done

#rm $tmp"/rand_seq.fasta"
#rm $tmp"/rand_unit.txt"
#rm $tmp"/rand_TRs.txt"

exit 0
