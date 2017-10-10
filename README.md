# mTR
Use Makefile to generate an executable file "mTR".

Execution
mTR -a|-r|-c <fasta file> (-pp option)

Arguments
-a: Perform All steps. Feed <fasta file>, detect repeats, and print <clustering results>.
  
-r: Find Repeats in reads only. Feed <fasta file> and print <reads with repeats>.
  
-c: Cluster reads with repeats. Feed <reads with repeats> and print <clustering results>.
  
-pp: Pretty print of reads with repeats for better readability and debugging
