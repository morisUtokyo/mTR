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
    void chaining(int print_alignment);
    void search_max(int print_alignment);
    void delete_set_of_alignments();
    extern void pretty_print_alignment(char *unit_string, int unit_len, int rep_start, int rep_end);
    extern void wrap_around_DP(
                               int *rep_unit, int unit_len,
                               int query_start, int query_end,
                               int *actual_start,   int *actual_end,
                               int *return_rep_len, int *return_freq_unit,
                               int *return_matches, int *return_mismatches,
                               int *return_insertions, int *return_deletions);
}
#endif
