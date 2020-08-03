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

// Key default parameters
#define MAX_INPUT_LENGTH 1000000 // The maximum length of each read
#define MIN_MATCH_RATIO 0.6      // The minimum threshold of match ratio　between the estimated repeat unit and the repeat in a given raw read
#define MIN_PERIOD 2            // Minimum period length
#define MAX_PERIOD 500          // Maximum period length
#define MIN_NUM_FREQ_UNIT 5     // The minimum threshold of number of units

#define ALIGNMENT_WIDTH_PRINTING 50
#define MAX_LEN_overlapping 10  // = MIN_PERIOD * MIN_NUM_FREQ_UNIT


// The following values are optimzed for a benchmark dataset.
#define MIN_WINDOW 5
#define MAX_WINDOW 10240

#define minKmer 5
#define maxKmer 15      // Increase this when no qualified repeats are found.
#define MAX_tiebreaks 1024
#define MIN_jaccard_index 0.98

#define BLK 4096        // Block size of input buffer.
#define WrapDPsize  200000000    // 200M  = repeat_unit_size (200) x length_of_repeats (1,000,000)

//  Choice of DeBruijn graph or progressive multiple alignment
#define ProgressiveMultipleAlignment 0
#define DeBruijnGraphSearch 1

#define count_maxKmer 6
#define PrimeMax 256019

//  Global variables in the heap
int Manhattan_Distance;         // The default setting is 1 in the main.
float min_match_ratio;
int *pow4;              // A table of pow(4, k)
int *min_coverage_for_missing_bases;
int *orgInputString;    // An input read string of length MAX_INPUT_LENGTH.
int *inputString;       // 4 decimal encoding of the input read string of length MAX_INPUT_LENGTH.
int *inputString_w_rand; // 4 decimal encoding of the input read string of length MAX_INPUT_LENGTH.
int *count;             // A table of size 4^k for counting sort.
int **freqNode;         // A hash table for storing the frequency of each node (node, frequency)
int *sortedString;      // Positions of 4-mers sorted wrt both 4-mers and their positions.

double *directional_index_tmp;  // For storing all DI values temporarily for a given w
double *directional_index;      // For storing a locally maximum DI
int *directional_index_end;     // The end position of the local maximum
int *directional_index_w;       // The window width w for the local maximum

int *vector0, *vector1, *vector2;
double *freq_interval_len;       // For computing frequency distribution of interval lengths
                       // The length is MAX_INPUT_LENGTH
int *WrapDP;            // 2D space for Wrap-around global alignment DP for handling tandem repeats
char *alignment_input;
char *alignment_symbols;
char *alignment_repeats;
// For printing the alignment of the predicted repeat unit with the input string
                        // The largest array, and the size is (MAX_PERIOD+1) * (MAX_INPUT_LENGTH+1)
int **consensus, **gaps;  // Space for consensus


typedef struct {        // MAX_ID_LENGTH + MAX_EPRIOD + 28*4 = 612 bytes
    int     ID;  // 0,1,2,...
    char    readID[BLK];
    int     inputLen;
    int     rep_start;
    int     rep_end;
    int     repeat_len;
    int     rep_period;
    int     Num_freq_unit;
    int     Num_matches;
    int     Num_mismatches;
    int     Num_insertions;
    int     Num_deletions;
    int     Kmer;
    int     match_gain;
    int     mismatch_penalty;
    int     indel_penalty;
    char    string[MAX_PERIOD];
    int     string_score[MAX_PERIOD];
    int     freq_2mer[16];
} repeat_in_read;

// External functions
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define DIFF(x, y) ((x) > (y) ? ((x) - (y)) : ((y) - (x)))

int handle_one_file(char *inputFile, int print_alignment);
void handle_one_read( char *readID, int inputLen, int read_cnt, int print_alignment);
void fill_directional_index_with_end(int DI_array_length, int inputLen, int random_string_length);
void init_inputString(int k, int query_start, int query_end, int inputLen);

int search_De_Bruijn_graph(int query_start, int query_end, repeat_in_read *rr);
void revise_representative_unit( repeat_in_read *rr);
void wrap_around_DP(int query_start, int query_end, repeat_in_read *rr);
void wrap_around_DP_sub( int query_start, int query_end, repeat_in_read *rr, int MATCH_GAIN, int MISMATCH_PENALTY, int INDEL_PENALTY );

void set_rr(repeat_in_read *rr_a, repeat_in_read *rr_b);
void clear_rr(repeat_in_read *rr_a);
void freq_2mer_array(int* val, int len, int *freq_2mer);
void print_one_repeat_in_read(repeat_in_read rr);
void print_freq(int rep_start, int rep_end, int rep_period, char* string, int inputLen, int k);

float time_all, time_memory, time_range, time_period, time_initialize_input_string, time_wrap_around_DP, time_count_table, time_chaining;
int query_counter;

//For debugging with #ifdef
//#define DEBUG_match_gain

