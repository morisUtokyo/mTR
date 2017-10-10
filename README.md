# mTR
Use Makefile to generate an executable file "mTR".

Execution
mTR -a|-r|-c <fasta file> (-pp option)

Arguments
-a: Perform All steps. Feed <fasta file>, detect repeats, and print <clustering results>.
-r: Find Repeats in reads only. Feed <fasta file> and print <reads with repeats>.
-c: Cluster reads with repeats. Feed <reads with repeats> and print <clustering results>.
-pp: Pretty print of reads with repeats for better readability and debugging


/*-----------------------------------------
Features associated with one read with repeats

repID         Identifier of a group with similar patterns of tandem repeats
ID             Identifier given by the program
readID        Identifier in the input fasta file
READ_length    Length of the focal read
REPEAT_start    Start position of a repeat estimated using Kadane's algorithm
end         End position of the estimated repeat
len(act)         Length of the actual repeat determined
UNIT len(est)   Length of the unit of the estimated repeat
len(act)           Length of the unit of the actual repeat
#Units        Number of unit occurrences in the actual repeat
Matches         Number of matches between the actual repeat and the tandem repeat
Mismacthes      Number of mismatches
Ins        Number of insertions
Del           Number of deletions
Kmer            Kmer used for calculating the actual repeat
Method        Method for computing the repeat unit: De Bruijn graph search or progressive multiple alignment
Unit_string     String of the unit of the actual repeat
2mer_freq    Frequency distribution of sixteen 2mers (AA,AC,AG, ...,TT)


Example
repID = 0,
ID = 0,
readID = 10_20_0_5_5_100_100_10.fasta_0,
READ_length = 397,
REPEAT_start = 100,
end = 293,
len(act) = 194,
UNIT_len(est) = 10,
len(act) = 10,
#Units = 19,
Matches = 0.964,
Mismacthes = 0.000,
Ins = 0.041,
Del = 0.036,
Kmer = 4,
Method = De Bruijn graph,
Unit_string= ACGTGCGGTA,
2mer_freq = 1,1,0,0,0,0,2,0,0,1,1,2,1,0,1,0,



