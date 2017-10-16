//
//  mTR.h
//  
//
//  Created by Shinichi Morishita on 2017/09/20.
//
//

// Key parameters

#define MAX_NUM_READS    1000000   // The maximum number of all reads in a given fasta file. Each entry need 1KB.
#define MAX_INPUT_LENGTH 200000  // The maximum length of each read
#define minKmer 3
#define maxKmer 11              // Increase this when no qualified repeats are found.
                            // Avoid using 15 =< k because ehe size is 4^maxKmer * 4B.
#define MIN_MATCH_RATIO 0.7      // The minimum threshold of match ratio
                            // between the estimated repeat unit and the repeat in a given raw read
#define BLK 1024                // Block size of input buffer.
#define MAX_PERIOD 500          // Maximum period length
#define MIN_REP_LEN_RATIO 0.6   // The threshold of the ratio of the actual repeat unit length to the estimated.
#define MIN_REP_LEN 20
#define MH_distance_threshold 0.2  // Two 2mer frequency distributions are identical if their Manhattan distance is less than or equal to this threshold.  A small threshold generates smaller groups of repeat units.

#define maxThresholdMultipleAlignment  20

#define MATCH_GAIN  1       // Parameters for global and wrap around alignment
#define MISMATCH_PENALTY -1
#define INDEL_PENALTY -1

//  Global variables stored in the heap
int *orgInputString;    // An input read string of length MAX_INPUT_LENGTH.
int *inputString;       // 4 decimal encoding of the input read string of length MAX_INPUT_LENGTH.
int *count;             // A table of size 4^k for counting sort.
int *sortedString;      // Positions of 4-mers sorted according to both 4-mers and their positions.
                        // The length is MAX_INPUT_LENGTH
float *freq_interval_len;       // For computing frequency distribution of interval lengths
                       // The length is MAX_INPUT_LENGTH
float *Kadane_val;      // List of values for Kadena's algorithm. The length is MAX_INPUT_LENGTH
int *max_starts;        // List of Starting positions of temporary maximums. The length is MAX_INPUT_LENGTH.
int *count_period_all;  // Frequency of individual periods of size 1 to MAX_PERIOD.
int *rep_unit_string;   // String of representative unit of the focal repeat. The length is MAX_PERIOD.
int *WrapDP;            // 2D space for Wrap-around global alignment DP for handling tandem repeats
                        // The largest array, and the size is (MAX_PERIOD+1) * (MAX_INPUT_LENGTH+1)

//  Arrays for storing results
#define MAX_ID_LENGTH 200
#define ProgressiveMultipleAlignment 0
#define DeBruijnGraphSearch 1

typedef struct {        // MAX_ID_LENGTH + MAX_EPRIOD + 28*4 = 612 bytes
    int     ID;  // 0,1,2,...
    char    readID[MAX_ID_LENGTH];
    int     inputLen;
    int     max_start;
    int     max_end;
    int     actual_repeat_len;
    int     rep_period;
    int     actual_rep_period;
    int     Num_freq_unit;
    int     Num_matches;
    int     Num_mismatches;
    int     Num_insertions;
    int     Num_deletions;
    int     Kmer;
    int     ConsensusMethod;     // 0 = progressive multiple alignment, 1 = De Bruijn graph search
    char    string[MAX_PERIOD];
    int     freq_2mer[16];
} repeat_in_read;

repeat_in_read *repeats_in_all_reads;

// For clustering repeats

#define MIN_NUM_repTR 1

typedef struct TR1{    // Sorted according to actual_rep_period, freq_2mer, and Num_freq_unit
    int     ID;     //  20 B
    int     actual_rep_period;
    int     freq_2mer[16];
    int     Num_freq_unit;
    int     repID;
    struct TR1 *repTR;      // This allows us to refer the struct itself.
    int     freq;
} TR;

TR      *TR_list;
TR   *repTR_list;

// External functions
#define MAX(a, b) ((a) > (b) ? (a) : (b))

int handle_one_file(char *inputFile);
void handle_one_read( char *readID, int inputLen, int read_cnt);

int progressive_multiple_alignment(
   int max_start, int max_end, int max_pos,
   int rep_period, int Kmer, int inputLen, int pow4k);
int search_De_Bruijn_graph(int max_pos, int rep_period, int inputLen, int pow4k_1);
void wrap_around_DP(int *rep_unit, int unit_len, int *rep, int rep_len,
   int *return_rep_len, int *return_freq_unit,
   int *return_matches, int *return_mismatches,
   int *return_insertions, int *return_deletions);

void freq_2mer_array(int* val, int len, int *freq_2mer);
void pretty_print_one_repeat_in_read(repeat_in_read rr);
void print_one_repeat_in_read(repeat_in_read rr);
int feed_rr_into_repeats_in_all_reads(char *inputFile);

void k_means_clustering(int read_cnt, int pretty_print);
void print_one_TR_with_read(TR aTR, int pretty_print);

// For debugging with #ifdef

//#define DEBUG_IO
//#define DEBUG_sorting
//#define DEBUG_algorithm
//#define DEBUG_algorithm_freq
//#define DEBUG_algorithm_Kadane
//#define DEBUG_algorithm_rep_freq
//#define DEBUG_algorithm_wrap_around_all
//#define DEBUG_algorithm_wrap_around
//#define DEBUG_progressive_multiple_alignment
//#define DEBUG_algorithm_print_rep_unit
//#define DEBUG_clustering1
//#define DEBUG_clustering2
//#define DEBUG_repeats_in_all_reads
//#define output_clustering
