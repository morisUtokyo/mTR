/*
Copyright (c) 2019, Shinichi Morishita
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
*/

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
                          char* string,
                          int*  string_score,
                          int   match_gain,
                          int   mismatch_penalty,
                          int   indel_penalty
                          );
    void chaining(int print_alignment);
    void search_max(int print_alignment);
    void delete_set_of_alignments();
    extern void pretty_print_alignment(char *unit_string, int unit_len, int rep_start, int rep_end, int MATCH_GAIN, int MISMATCH_PENALTY, int INDEL_PENALTY);
    extern void print_freq(int rep_start, int rep_end, int rep_period, char* string, int inputLen, int k);
}
#endif
