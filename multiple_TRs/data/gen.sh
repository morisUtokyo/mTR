#!/bin/bash

for i in 3 5 10 20 50; do
for j in 3 5 10 20 50; do
if [ $i -lt $j ]; then
echo $i"_"$j
../rand_multi_seq $i"_"$j"_set.txt" $i"_"$j".fasta" $i"_"$j"_unit.txt"
fi
done
done

exit 0
