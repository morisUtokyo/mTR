# mTR
Use Makefile to generate an executable file "mTR".

Execution
mTR -m|-r <fasta file> 

Arguments

-r: Find one tandem repeat in each read. 
-m: Output multiple occurrences of different tandem repeats.

/*-----------------------------------------
Features associated with one read with repeats

read identifier        
read length   
repeat start position
repeat end position 
repeat length      
unit length    
number of unit occurrences 
Matches        Number of matches between the actual repeat and the tandem repeat
Match ratio   Ratio of the number of matches to the length of actual repeat determined
Mismacthes  Number of mismatches
Insertions      Number of insertions
Deletions       Number of deletions
Kmer             Kmer used for calculating the actual repeat
Method         Method for computing the repeat unit: De Bruijn graph search (1) or progressive multiple alignment (0)
unit string


Example
read identifier = 10_20_0_5_5_100_100_10.fasta_0,
read length = 1200,
repeat start position = 28,
repeat end position  = 1090,
repeat length = 1063,
unit length = 101,
number of unit occurrences = 10,
Matches  = 852,
Match ratio = 0.8015
Mismacthes = 137,
Insertions  = 101,
Deletions  = 74,
Kmer  = 7,
Method   = De Bruijn graph (1) or Progressive multiple alignment (0)
unit string = GTGACCCCGCGGTAGCATAGCGGCTGAATGCCAGATACACATGTCAAGTCGTGGTGGCCCCGGTACTGGCTTCCGAGGCTGTAGAAAACCTATTCCCAGCG,



