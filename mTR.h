//
//  mTR.h
//  
//
//  Created by Shinichi Morishita
//
//

// Key parameters

#define MAX_INPUT_LENGTH 200000  // The maximum length of each read
#define MIN_MATCH_RATIO 0.7      // The minimum threshold of match ratio　between the estimated repeat unit and the repeat in a given raw read

// The following values are optimzed for a benchmark dataset.
#define minKmer 5
#define maxKmer1 9
#define maxKmer2 11              // Increase this when no qualified repeats are found.
#define BLK 1024                // Block size of input buffer.
#define MAX_PERIOD 500          // Maximum period length
#define MIN_PERIOD 2            // Minimum period length
#define MIN_NUM_FREQ_UNIT 5     // The minimum threshold of number of units
#define WrapDPsize  20000000    // 10M  > repeat_unit_size (100) x num_of_repeats (100,000)
// Parameters for global and wrap around alignment
#define MATCH_GAIN  1
#define MISMATCH_PENALTY -1
#define INDEL_PENALTY -1
//  Choice of DeBruijn graph or progressive multiple alignment
#define ProgressiveMultipleAlignment 0
#define DeBruijnGraphSearch 1

//  Global variables in the heap
int *pow4;              // A table of pow(4, k)
int *orgInputString;    // An input read string of length MAX_INPUT_LENGTH.
int *inputString;       // 4 decimal encoding of the input read string of length MAX_INPUT_LENGTH.
int *count;             // A table of size 4^k for counting sort.
int *sortedString;      // Positions of 4-mers sorted wrt both 4-mers and their positions.
double *directional_index;  // For detecting boundaries of tandem repeats

int *vector0, *vector1, *vector2;
double *freq_interval_len;       // For computing frequency distribution of interval lengths
                       // The length is MAX_INPUT_LENGTH

int *count_period_all;  // Frequency of individual periods of size 1 to MAX_PERIOD.
int *rep_unit_string;   // String of representative unit of the focal repeat. The length is MAX_PERIOD.
int *WrapDP;            // 2D space for Wrap-around global alignment DP for handling tandem repeats
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
    int     ConsensusMethod;     // 0 = progressive multiple alignment, 1 = De Bruijn graph search
    char    string[MAX_PERIOD];
    int     predicted_rep_period;
    int     freq_2mer[16];
} repeat_in_read;

repeat_in_read *RRs;


// External functions
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

int handle_one_file(char *inputFile, int print_multiple_TR);
void handle_one_read( char *readID, int inputLen, int read_cnt, int print_multiple_TR);

int progressive_multiple_alignment(
   int max_start, int max_end, int max_pos,
   int rep_period, int Kmer, int inputLen, int pow4k);
int search_De_Bruijn_graph(int max_start, int max_end, int max_pos,
   int inputLen, int pow4k_1);
void wrap_around_DP(
   int *rep_unit, int unit_len, int *rep, int rep_len,
   int *actual_start,   int *actual_end,
   int *return_rep_len, int *return_freq_unit,
   int *return_matches, int *return_mismatches,
   int *return_insertions, int *return_deletions);

void freq_2mer_array(int* val, int len, int *freq_2mer);
void print_one_repeat_in_read(repeat_in_read rr);

// For debugging with #ifdef

//#define DEBUG_algorithm_wrap_around_all
//#define DEBUG_algorithm_wrap_around
//#define DEBUG_progressive_multiple_alignment
//#define DEBUG_incremental

//#define DEBUG_window_kmer

float time_all, time_memory, time_range, time_period, time_predicted_rep_period_and_max_position, time_wrap_around_DP, time_search_De_Bruijn_graph, time_progressive_multiple_alignment, time_count_table;

