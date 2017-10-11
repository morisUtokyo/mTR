//
//  handle_one_read.c
//  
//
//  Created by Shinichi Morishita on 2017/10/06.
//
// vc++ disable 4996
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"


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

void print_4_decimal_array(int* val, int len, char *return_string){
    strcpy(return_string, "");
    for(int i=0; i<len; i++){
        switch(val[i]){
            case 0: strcat(return_string, "A"); break;
            case 1: strcat(return_string, "C"); break;
            case 2: strcat(return_string, "G"); break;
            case 3: strcat(return_string, "T"); break;
            default: fprintf(stderr, "fatal error: input char %i", val[i]); exit(EXIT_FAILURE);
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


void handle_one_read_with_a_Kmer(
                                 char *readID,
                                 int inputLen,
                                 int Kmer,
                                 repeat_in_read *rr
                                 ){
    int i;
    
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
    for(i = 0; i < pow4k; i++){ count[i] = 0;}  // Initialization
    for(i = 0; i < inputLen; i++){ count[ inputString[i] ]++; }    // Perform counting
    for(i = 1; i < pow4k; i++){ count[i] = count[i-1] + count[i]; }
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
    //  Compute the frequency distribution of interval lengths by counting sort
    //---------------------------------------------------------------------------
    
    for(int j=0; j<MAX_INPUT_LENGTH; j++){  // Initialization
        freq_interval_len[j] = 0;
    }
    for(int j=0; j<inputLen; j++){
        start = count[ inputString[j] ];
        if(inputString[j] == pow4k-1){
            end = inputLen;
        }else{
            end = count[ inputString[j]+1 ];
        }
        freq_interval_len[j] = (float)(end - start);
    }
    
    
    // Compute the actual average and standard deviation
    float avg = 0;
    float sd = 0;
    for(int j=0; j <inputLen; j++){
        avg += freq_interval_len[j];
    }
    avg = avg/inputLen;
    for(int j=0; j <inputLen; j++){
        sd += (freq_interval_len[j] - avg)*(freq_interval_len[j] - avg);
    }
    sd = sqrtf(sd/inputLen);
    
    // Assuming the random distribution of k-mers, computethe averaga and standard deviation
    float p = (float)1/pow4k;
    float random_avg = inputLen * p;
    float random_sd  = sqrtf(inputLen * p * (1-p));
    
    for(int j=0; j <inputLen; j++){
        freq_interval_len[j] -= (1 + random_avg + 3 * random_sd);
                // 1 means the count of itself.  4 * standard deviation is used.
        // freq_interval_len[j] = freq_interval_len[j] - avg;
                // This ad hoc approach was avoided.
    }
    
#ifdef DEBUG_algorithm_freq
    printf("avg=%f\t sd=%f p=%f pow4k=%i len=%i random_avg=%f random_sd=%f\n", avg, sd, p, pow4k, inputLen, random_avg, random_sd);
    //printf("position\tfreq\n");
    //for(int j=0; j<inputLen; j++){
    //    printf("%i\t%i\n", j, (int)freq_interval_len[j]);
    //}
    //printf("\n");
#endif
    
    //---------------------------------------------------------------------------
    //  Determine a highly repetitive region as an optimal range using Kadane's algorithm
    //---------------------------------------------------------------------------
    
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
    int max_end = 0;
    for(int j=1; j<inputLen; j++){
        if(max_val < Kadane_val[j]){
            max_val = Kadane_val[j];
            max_end = j;
        }
    }
    int max_start = max_starts[max_end];
    
#ifdef DEBUG_algorithm_Kadane
    printf("position\tfreq\n");
    for(int j=max_start; j<max_end; j++){
        printf("%i\t%i\n", j, (int)freq_interval_len[j]);
    }
    printf("\n");
    printf("avg=%f\t sd=%f\n", avg, sd);
#endif
    
    //---------------------------------------------------------------------------
    // Compute the periodicity with the maximum count
    //---------------------------------------------------------------------------
    
    // Initialization
    for(i = 0; i < MAX_PERIOD; i++){
        count_period_all[i] = 0;
    }
    for(int i=max_start; i <= max_end; i++){
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
    // Compute the representative unit string by using a progressive multiple alignment or
    // by traversing the De Bruijn graph of all k-mers in a greedy manner
    //---------------------------------------------------------------------------
    
    // Initialization
    int actual_rep_period;
    for(i = 0; i < MAX_PERIOD; i++){
        rep_unit_string[i] = 0;
    }
    //　Locate the initial position (k-mer) with the maximum frequency
    int max_pos = max_start;
    max_freq   = (int)freq_interval_len[max_pos];
    for(int j=max_start; j <= max_end; j++){
        if(max_freq < freq_interval_len[j]){
            max_pos  = j;
            max_freq = (int)freq_interval_len[max_pos];
        }
    }
    
    // First, traverse the De Bruijn graph of all k-mers in a greedy manner
    actual_rep_period =
        search_De_Bruijn_graph(max_pos, rep_period, inputLen, pow4k_1);

    //  When a repeat unit is found, actual_rep_period > 0, and = 0 otherwise.
    int ConsensusMethod;
    if(actual_rep_period > 0){
        // A De Bruijin graph search is successful.
        ConsensusMethod = DeBruijnGraphSearch;
    }else{
        // If the De Bruijn graph search fails, try a progressive  multiple alignment.
        ConsensusMethod = ProgressiveMultipleAlignment;
        actual_rep_period = progressive_multiple_alignment(
          max_start, max_end, max_pos, rep_period, Kmer, inputLen, pow4k);
    }

    
    //---------------------------------------------------------------------------
    // Compute the accuracy of the representative unit string by wrap-around DP
    //---------------------------------------------------------------------------
    
    strcpy( rr->readID, readID);
    rr->inputLen  = inputLen;
    rr->max_start = max_start;
    rr->max_end   = max_end;
    rr->rep_period= rep_period;
    rr->actual_rep_period = actual_rep_period;
    rr->Kmer      = Kmer;
    rr->Num_freq_unit     = 0;
    rr->ConsensusMethod = ConsensusMethod;
    
    if(actual_rep_period > 0){
        int actual_repeat_len, Num_freq_unit, Num_matches, Num_mismatches, Num_insertions, Num_deletions;
        
        wrap_around_DP(rep_unit_string,
                       actual_rep_period,
                       &orgInputString[max_start],
                       (max_end - max_start + 1),
                       &actual_repeat_len,  &Num_freq_unit,
                       &Num_matches,        &Num_mismatches,
                       &Num_insertions,     &Num_deletions);
        
        rr->Num_freq_unit     = Num_freq_unit;
        rr->Num_matches       = Num_matches;
        rr->Num_mismatches    = Num_mismatches;
        rr->Num_insertions    = Num_insertions;
        rr->Num_deletions     = Num_deletions;
        
        rr->actual_repeat_len = actual_repeat_len;
        print_4_decimal_array(rep_unit_string, actual_rep_period, rr->string);
        freq_2mer_array(rep_unit_string, actual_rep_period, rr->freq_2mer);
    }
}

void clear_rr(repeat_in_read *rr_a){
    rr_a->ID                = -1;
    strcpy( rr_a->readID, "");
    rr_a->inputLen          = -1;
    rr_a->max_start         = -1;
    rr_a->max_end           = -1;
    rr_a->actual_repeat_len = -1;
    rr_a->rep_period        = -1;
    rr_a->actual_rep_period = -1;
    strcpy( rr_a->string, "");
    for(int i=0; i<16; i++){ rr_a->freq_2mer[i] = -1; }
    rr_a->Num_freq_unit     = -1;
    rr_a->Num_matches       = -1;
    rr_a->Num_mismatches    = -1;
    rr_a->Num_insertions    = -1;
    rr_a->Num_deletions     = -1;
    rr_a->Kmer              = -1;
    rr_a->ConsensusMethod     = -1;
}

void cpy_rr(repeat_in_read *rr_a, repeat_in_read *rr_b){
    rr_a->ID                = rr_b->ID;
    strcpy( rr_a->readID, rr_b->readID);
    rr_a->inputLen          = rr_b->inputLen;
    rr_a->max_start         = rr_b->max_start;
    rr_a->max_end           = rr_b->max_end;
    rr_a->actual_repeat_len = rr_b->actual_repeat_len;
    rr_a->rep_period        = rr_b->rep_period;
    rr_a->actual_rep_period = rr_b->actual_rep_period;
    strcpy( rr_a->string, rr_b->string);
    for(int i=0; i<16; i++){ rr_a->freq_2mer[i] = rr_b->freq_2mer[i]; }
    rr_a->Num_freq_unit     = rr_b->Num_freq_unit;
    rr_a->Num_matches       = rr_b->Num_matches;
    rr_a->Num_mismatches    = rr_b->Num_mismatches;
    rr_a->Num_insertions    = rr_b->Num_insertions;
    rr_a->Num_deletions     = rr_b->Num_deletions;
    rr_a->Kmer              = rr_b->Kmer;
    rr_a->ConsensusMethod   = rr_b->ConsensusMethod;
}


void handle_one_read(char *readID, int inputLen, int read_cnt){
    
    // Put into repeats_in_all_reads any repeat instance that may not meet the conditions on MIN_REP_LEN and MIN_MATCH_RATIO. Remove unqualified instances later.
    int max_matches = -1;  // Initial value
    char message[1000] = "";
    repeat_in_read max_rr, rr;
    
    for(int k = minKmer; k <= maxKmer; k++){
        clear_rr(&rr);
        handle_one_read_with_a_Kmer(readID, inputLen, k, &rr);
        
        if( max_matches < rr.Num_matches ){
            max_matches = rr.Num_matches;
            cpy_rr(&max_rr, &rr);
            max_rr.ID = read_cnt;
        }
    }
    cpy_rr(&repeats_in_all_reads[max_rr.ID], &max_rr);
}

