//
//  handle_one_read.c
//  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

// Probability of x in normal distribution
float p_nd(float x, float avg, float sd){
    float p = (float)exp( -pow(x-avg,2) / (2 * pow(sd,2) ) ) / sqrt( 2 * M_PI * pow(sd,2) );
    return(p);
}

void print_4_decimals(int val, int len){
    if(len > 0){
        print_4_decimals(val/4, len-1);
        switch(val%4){
            case 0: printf("A"); break;
            case 1: printf("C"); break;
            case 2: printf("G"); break;
            case 3: printf("T"); break;
            default: fprintf(stderr, "fatal input char %i", val%4); exit(EXIT_FAILURE);
        }
    }
}

void freq_2mer_array(int* val, int len, int *freq_2mer){
    for(int i=0; i<16; i++){
        freq_2mer[i] = 0;
    }
    for(int i=1; i<len; i++){
        freq_2mer[val[i-1]*4+val[i]]++;
    }
    // wrap around and concatenate the last and first characters
    freq_2mer[val[len-1]*4+val[0]]++;
}

// Return 1 if max_rr < rr, 0 otherwise.
int better_rr(repeat_in_read *max_rr, repeat_in_read *rr){
    
    float max_rr_match_ratio = 0;
    if(max_rr->actual_repeat_len > 0){
        max_rr_match_ratio = (float)max_rr->Num_matches / max_rr->actual_repeat_len;
    }
    float rr_match_ratio = 0;
    if(rr->actual_repeat_len > 0){
        rr_match_ratio = (float)rr->Num_matches / rr->actual_repeat_len;
    }
    
    if( MIN_MATCH_RATIO < rr_match_ratio &&
       max_rr_match_ratio < rr_match_ratio &&
       max_rr->Num_matches < rr->Num_matches &&
       MIN_REP_LEN < (rr->actual_rep_period * rr->Num_freq_unit) &&
       1 < rr->Num_freq_unit )
    {
        return(1);
    }else{
        return(0);
    }
}

void substitute_rr(repeat_in_read *rr_a, repeat_in_read *rr_b){
    rr_a->ID                = rr_b->ID;
    strcpy( rr_a->readID, rr_b->readID);
    rr_a->inputLen          = rr_b->inputLen;
    rr_a->max_start         = rr_b->max_start;
    rr_a->max_end           = rr_b->max_end;
    rr_a->actual_repeat_len = rr_b->actual_repeat_len;
    rr_a->rep_period        = rr_b->rep_period;
    rr_a->actual_rep_period = rr_b->actual_rep_period;
    rr_a->Num_freq_unit     = rr_b->Num_freq_unit;
    rr_a->Num_matches       = rr_b->Num_matches;
    rr_a->Num_mismatches    = rr_b->Num_mismatches;
    rr_a->Num_insertions    = rr_b->Num_insertions;
    rr_a->Num_deletions     = rr_b->Num_deletions;
    rr_a->Kmer              = rr_b->Kmer;
    rr_a->ConsensusMethod   = rr_b->ConsensusMethod;
    strcpy( rr_a->string, rr_b->string);
    for(int i=0; i<16; i++){ rr_a->freq_2mer[i] = rr_b->freq_2mer[i]; }
}

repeat_in_read handle_one_read_with_a_Kmer( char *readID, int inputLen, int Kmer ){
    int pow4k_1 = 1;    // 4^{k-1}  e.g., 4^(4-1) = 64
    for(int i=0; i<(Kmer-1); i++){ pow4k_1 = 4 * pow4k_1; }
    int pow4k = 4 * pow4k_1;  // 4 * 4^{k-1} = 4^k
    
    // Encode the raw input string into 4 decimals of length k
    for(int i=0; i<inputLen; i++){
        inputString[i] = orgInputString[i];
    }
    int tmp = 0;
    for(int i=0; i<(Kmer-1); i++){
        tmp = 4 * tmp + inputString[i];  // compute 4 decimal of length k-1
    }
    for(int i=0; i<(inputLen-Kmer+1); i++){
        inputString[i] = 4 * tmp + inputString[i+Kmer-1];
        tmp = inputString[i] % pow4k_1; //　remainder, compute 4 decimal of length k-1
        if(tmp < 0){
            fprintf(stderr, "fatal error at %i\t %i \t%i\n", i, inputString[i], pow4k_1);
            exit(EXIT_FAILURE);
        }
    }
    inputLen = inputLen - Kmer + 1;
    
#ifdef DEBUG_IO
    printf("Input length of k decimals = %i\n", inputLen);
    for(int i=0; i<inputLen; i++){
        printf("%i ", inputString[i]);
    }
    printf("\n");
#endif
    
    //---------------------------------------------------------------------------
    // Sort k-mers using counting sort and store the positions of k-mers in sortedString
    //---------------------------------------------------------------------------
    for(int i = 0; i < pow4k; i++){ count[i] = 0;}  // Initialization
    for(int i = 0; i < inputLen; i++){ count[ inputString[i] ]++; }    // Perform counting
    for(int i = 1; i < pow4k; i++){ count[i] = count[i-1] + count[i]; }
    for(int i=0; i<MAX_INPUT_LENGTH; i++){ sortedString[i] = 0; } // Initialization
    for(int i=inputLen-1; 0 <= i; i--){
        sortedString[ --count[inputString[i]] ] = i;
    }
    
#ifdef DEBUG_sorting
    printf("Sorted \n");
    for(int i=0; i<inputLen; i++){
        printf("%i ", sortedString[i]);
    }
    printf("\n");
    
    printf("4dec\t counts\n");
    for(i=1; i<pow4k; i++){
        printf("%i\t %i\n", i, count[i] );
    }
#endif
    
    int start, end;
    
#ifdef DEBUG_algorithm
    printf("Display M(i,j) \n");
    int* oneRow = (int *)malloc( sizeof(int) * inputLen);
    for(int i=0; i<inputLen; i++){
        for(int j=0; j<inputLen; j++){
            oneRow[j]=0;
        }  // Initialize one row
        start = count[ inputString[i] ];
        if(inputString[i] == pow4k-1){
            end = inputLen;
        }else{
            end   = count[ inputString[i]+1 ];
        }
        for(int j = start; j < end; j++){
            oneRow[ sortedString[j] ] = 1;
        }
        for(int j=0; j<inputLen; j++){
            printf("%i\t", oneRow[j]);
        }
        printf("\n");
    }
    free(oneRow);
#endif
    
    //---------------------------------------------------------------------------
    // Compute the periodicity with the maximum count in the input string
    //---------------------------------------------------------------------------
    
    // Initialization
    for(int i = 0; i < MAX_PERIOD; i++){
        count_period_all[i] = 0;
    }
    for(int i=0; i < inputLen; i++){
        start = count[inputString[i]];
        if(inputString[i] == pow4k-1){
            end = inputLen;
        }else{
            end = count[ inputString[i]+1 ];
        }
        for(int j = start; j < (end-1); j++){
            int aPeriod = sortedString[j+1] - sortedString[j];
            if(aPeriod < MAX_PERIOD){
                count_period_all[ aPeriod ]++;
            }
        }
    }
    int rep_period = 2;
    int max_freq = 0;
    for(int p = 2; p < MAX_PERIOD; p++){
        if(max_freq < count_period_all[p]){
            rep_period = p;
            max_freq = count_period_all[p];
        }
    }
    
    //---------------------------------------------------------------------------
    //  Compute the frequency distribution of interval lengths by counting sort
    //---------------------------------------------------------------------------
    
    for(int j=0; j<MAX_INPUT_LENGTH; j++){  // Initialization
        freq_interval_len[j] = 0;
    }
    
    int band_width = rep_period * 10;  // 10 performed best among {3,5,10,20,50}
    // To make sure the independence of the input length, examine the band within 10 times the estimated repeat unit length from the diagonal
    
    for(int j=0; j<inputLen; j++){
        start = count[ inputString[j] ];
        if(inputString[j] == pow4k-1){
            end = inputLen;
        }else{
            end = count[ inputString[j]+1 ];
        }

        int cnt = 0;
        for(int h = start; h < end; h++){
            if( (j - band_width) <= sortedString[h] && sortedString[h] <= (j + band_width)){
                cnt++;
            }
        }
        freq_interval_len[j] = (float)(cnt - 1);
        // freq_interval_len[j] = (float)(end - start) - 1;
        // Remove the k-mer at the j-th position and reduce the frequency by 1
    }
    
    // Assuming the random distribution of k-mers, computethe averaga and standard deviation
    float p = (float)1/pow4k;
    float random_avg = 2 * band_width * p;
    float random_sd  = sqrtf(2 * band_width * p * (1-p));
    
    for(int j=0; j <inputLen; j++){
        /*
         // Probablistic model did not perform well
        double p_value = p_nd( freq_interval_len[j], random_avg, random_sd);
        double p_threshold = 0.05;
        
        if(p_value < (p_threshold / inputLen) ){  // Bonferroni correction as we examine "inputLen" positions.
        //if(p_value < p_threshold  ){
            freq_interval_len[j] = Kmer;
        }else{
            freq_interval_len[j] = -1;
        }
         */
        
        freq_interval_len[j] -= (random_avg + 1 * random_sd);
                    // best prediction among {0,1,2,3}
    }

    //---------------------------------------------------------------------------
    //  Determine a highly repetitive region as an optimal range using Kadane's algorithm
    //---------------------------------------------------------------------------

    int max_start, max_end;

    // Initialization
    for(int j=0; j<MAX_INPUT_LENGTH; j++){
        Kadane_val[j] = 0;
        max_starts[j] = 0;
    }
    for(int j=0; j<inputLen; j++){
        if(j == 1){
            Kadane_val[j] = freq_interval_len[j];
            max_starts[j] = j;
        }else{
            if(freq_interval_len[j] < Kadane_val[j-1]+freq_interval_len[j]){
                Kadane_val[j] = Kadane_val[j-1]+freq_interval_len[j];
                max_starts[j] = max_starts[j-1];
            }else{
                Kadane_val[j] = freq_interval_len[j];
                max_starts[j] = j;
            }
        }
    }
    // backtracing
    float max_val = Kadane_val[0];

    max_end = 0;
    for(int j=1; j<inputLen; j++){
        if(max_val < Kadane_val[j]){
            max_val = Kadane_val[j];
            max_end = j;
        }
    }
    max_start = max_starts[max_end];
    

#ifdef DEBUG_algorithm_Kadane
    printf("position\tfreq\n");
    for(int j=max_start; j<max_end; j++){
        printf("%i\t%i\n", j, (int)freq_interval_len[j]);
    }
    printf("\n");
#endif

    
    //---------------------------------------------------------------------------
    // Compute the representative unit string by using a progressive multiple alignment or
    // by traversing the De Bruijn graph of all k-mers in a greedy manner
    //---------------------------------------------------------------------------

    repeat_in_read tmp_rr;
    strcpy( tmp_rr.readID, readID);
    tmp_rr.inputLen  = inputLen;
    tmp_rr.max_start = max_start;
    tmp_rr.max_end   = max_end;
    tmp_rr.rep_period= rep_period;
    tmp_rr.Kmer      = Kmer;
    tmp_rr.actual_repeat_len = 0;
    tmp_rr.Num_freq_unit     = 0;

    repeat_in_read max_rr, rr;
    substitute_rr(&rr, &tmp_rr);
    
    // Set max_rr to the NULL value
    max_rr.Num_matches = -1;
    max_rr.actual_repeat_len = 0;
    
    // First, traverse the De Bruijn graph of all k-mers in a greedy manner
    search_De_Bruijn_graph(pow4k, &tmp_rr, &rr);
    if( better_rr(&max_rr, &rr) == 1 ){
        substitute_rr(&max_rr, &rr);
    }else{
#ifdef USE_an_additional_progressive_multiple_alignment
        // If De Bruijn search fails, use a progressive multiple alignment
        progressive_multiple_alignment(pow4k, &tmp_rr, &rr);
        if( better_rr(&max_rr, &rr) == 1){
            substitute_rr(&max_rr, &rr);
        }
#endif
    }
    return(max_rr);
}

void handle_one_read(char *readID, int inputLen, int read_cnt){
    
    // Put into repeats_in_all_reads any repeat instance that may not meet the conditions on MIN_REP_LEN and MIN_MATCH_RATIO. Remove unqualified instances later.
    
    repeat_in_read max_rr, rr;
    
    // Set max_rr to the NULL value
    max_rr.Num_matches = -1;
    max_rr.actual_repeat_len = 0;
    
    for(int k = minKmer; k <= maxKmer; k++){
        rr = handle_one_read_with_a_Kmer(readID, inputLen, k);
        //print_one_repeat_in_read(rr);
        if( better_rr(&max_rr, &rr) == 1){
            substitute_rr(&max_rr, &rr);   // max_rr = rr;
            max_rr.ID = read_cnt;
        }
        // rep_period is more accurate than actual_rep_period is for some kmer instances.
        if(     (max_rr.rep_period < 8   && k == maxKmer_lt_8)
           ||   (max_rr.rep_period < 20  && k == maxKmer_lt_20)
           ||   (max_rr.rep_period < 70  && k == maxKmer_lt_70)
           ){
            break;
        }
    }
    substitute_rr( &repeats_in_all_reads[max_rr.ID], &max_rr);
}

