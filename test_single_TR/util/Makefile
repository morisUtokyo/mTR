# Makefile
CXXFLAGS = -std=c++11
all: rand_seq count_match comp_mTR_DP
comp_mTR_DP: comp_mTR_DP.cpp
    c++ -std=c++11 -o comp_mTR_DP -Wcpp comp_mTR_DP.cpp
count_match: count_match.cpp
    c++ -std=c++11 -o count_match -Wcpp count_match.cpp
rand_seq: rand_seq.cpp
    c++ -std=c++11 -o rand_seq -Wcpp rand_seq.cpp
clean:
	rm -f rand_seq count_match comp_mTR_DP
