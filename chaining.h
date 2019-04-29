#ifdef __cplusplus
extern "C" {
    void construct_set_of_alignments();
    void insert_an_alignment_into_set(char* readID,
                          int   inputLen,
                          int   rep_start,
                          int   rep_end,
                          int   repeat_len,
                          int   rep_period,
                          int   Num_freq_unit,
                          int   Num_matches,
                          int   Num_mismatches,
                          int   Num_insertions,
                          int   Num_deletions,
                          int   Kmer,
                          int   ConsensusMethod,
                          char* string );
    void chaining();
    void search_max();
    void delete_set_of_alignments();
}
#endif
