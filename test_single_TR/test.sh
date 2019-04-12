#!/bin/bash

# executable modules
mTR=../mTR
rand_seq=util/rand_seq
count_match=util/count_match
# temporary directory
tmp=tmp
# error pattern
substitution_rate=10
insertion_rate=3
deletion_rate=8
# number of sequences
num_rand_seqs=1000
# Length of random sequences before and after the tandem repeat
rand_seq_len=1000


echo "This program generates a fasta file with "$num_rand_seqs" sequences with tandem repeats according to the following error pattern, predicts the unit string for each sequence, and outputs the number of correct answers in the prediction."

echo "substitution rate = " $substitution_rate "(%)"
echo "insertion rate    = " $insertion_rate "(%)"
echo "deletion rate     = " $deletion_rate "(%)"

# j frequency of units
for j in 10 30 50 100 200; do
# i unit length
for i in 2 3 4 5 6 7 8 9 10 50 100 200; do
echo "------------------"
echo "Frequency of units = "$j", Length of units = "$i
$rand_seq $tmp"/rand_seq.fasta" $tmp"/rand_unit.txt" $i $j $substitution_rate $insertion_rate $deletion_rate $rand_seq_len $rand_seq_len $num_rand_seqs
$mTR -s $tmp"/rand_seq.fasta" > $tmp"/rand_TRs.txt"
echo "Number of correct answers for" $num_rand_seqs "sequences with tandem repeats"
$count_match $tmp"/rand_TRs.txt" $tmp"/rand_unit.txt"
done
done

rm $tmp"/rand_seq.fasta"
rm $tmp"/rand_unit.txt"
rm $tmp"/rand_TRs.txt"

exit 0
