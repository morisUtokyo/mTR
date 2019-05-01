//
//  handle_one_read.c
//  
//
//  Created by Shinichi Morishita
//

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"
#include "chaining.h"

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

char dec2char(int val){
    char return_char;
    switch(val){
        case 0: return_char = 'A'; break;
        case 1: return_char = 'C'; break;
        case 2: return_char = 'G'; break;
        case 3: return_char = 'T'; break;
        default: fprintf(stderr, "fatal input char %i", val); exit(EXIT_FAILURE);
    }
    return(return_char);
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




void find_tandem_repeat_sub(int query_start, int query_end, char *readID, int inputLen, int k, repeat_in_read *rr ){
    
    //---------------------------------------------------------------------------
    // In the range [query_start, query_end],
    // compute the representative unit string
    // by traversing the De Bruijn graph of all k-mers in a greedy manner
    //---------------------------------------------------------------------------
    
    //int actual_rep_period = search_De_Bruijn_graph(query_start, query_end, k, inputLen, pow4[k-1] );
    int actual_rep_period = search_De_Bruijn_graph(query_start, query_end, inputLen, k);
    
    if(actual_rep_period < MIN_PERIOD){
        return;
    }
    
    int ConsensusMethod = DeBruijnGraphSearch;
    /*
     //  If a repeat unit is found, actual_rep_period > 0, and = 0 otherwise.
     if(actual_rep_period > 0){
     // A De Bruijin graph search is successful.
     ConsensusMethod = DeBruijnGraphSearch;
     }else{
     // If the De Bruijn graph search fails, try a progressive  multiple alignment.
     ConsensusMethod = ProgressiveMultipleAlignment;
     
     gettimeofday(&s, NULL);
     actual_rep_period = progressive_multiple_alignment( query_start, query_end, max_pos, predicted_rep_period, k, inputLen, pow4[k] );
     gettimeofday(&e, NULL);
     time_progressive_multiple_alignment
     += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
     }
     */
    
    //---------------------------------------------------------------------------
    // Compute the accuracy of the representative unit string by wrap-around DP
    //---------------------------------------------------------------------------
    
    struct timeval s, e;
    
    if(actual_rep_period * (query_end - query_start + 1) > WrapDPsize){
        fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
        clear_rr(rr);
    }else if(actual_rep_period > 0){
        int actual_start, actual_end, actual_repeat_len, Num_freq_unit, Num_matches, Num_mismatches, Num_insertions, Num_deletions;
        
        
        gettimeofday(&s, NULL);
        wrap_around_DP(rep_unit_string,
                       actual_rep_period,
                       &orgInputString[query_start],
                       (query_end - query_start + 1),
                       &actual_start,       &actual_end,
                       &actual_repeat_len,  &Num_freq_unit,
                       &Num_matches,        &Num_mismatches,
                       &Num_insertions,     &Num_deletions);
        gettimeofday(&e, NULL);
        time_wrap_around_DP += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
        
        strcpy( rr->readID, readID);
        rr->inputLen           = inputLen;
        rr->rep_start          = query_start + actual_start;
        rr->rep_end            = query_start + actual_end;
        rr->repeat_len         = actual_repeat_len;
        //rr->predicted_rep_period = predicted_rep_period;
        rr->rep_period         = actual_rep_period;
        rr->Kmer               = k;
        rr->Num_freq_unit      = 0;
        rr->ConsensusMethod    = ConsensusMethod;
        rr->Num_freq_unit      = Num_freq_unit;
        rr->Num_matches        = Num_matches;
        rr->Num_mismatches     = Num_mismatches;
        rr->Num_insertions     = Num_insertions;
        rr->Num_deletions      = Num_deletions;
        
        print_4_decimal_array(rep_unit_string, actual_rep_period, rr->string);
        freq_2mer_array(rep_unit_string, actual_rep_period, rr->freq_2mer);
    }
}

void init_inputString(int k, int inputLen){
    
    init_genrand(0);
    
    for(int i=0; i<inputLen; i++){ inputString[i] = (int)(genrand_int32()%4); }
    //for(int i=0; i<MAX_INPUT_LENGTH; i++){ inputString[i] = (int)(genrand_int32()%4); }
    
    for(int i=0; i<inputLen; i++){
        inputString[i] = orgInputString[i];
    }
    int tmp = 0;
    for(int i=0; i<(k-1); i++){
        tmp = 4 * tmp + inputString[i];  // compute 4 decimal of first k-1 letters
    }
    for(int i=0; i<(inputLen-k+1); i++){
        inputString[i] = 4 * tmp + inputString[i+k-1];
        tmp = inputString[i] % pow4[k-1]; //　remainder, compute 4 decimal of length k-1
        if(tmp < 0){
            fprintf(stderr, "fatal error at %i\t %i \t%i\n", i, inputString[i], pow4[k-1] );
            exit(EXIT_FAILURE);
        }
    }
}

void find_tandem_repeat(int query_start, int query_end, int w, char *readID, int inputLen, repeat_in_read *rr, repeat_in_read *tmp_rr ){
    
    // A heuristic rule for setting the k-mer range for a de Bruijn graph to an nearly optimal one
    int min_k, max_k;
    if(w < 100){                // 3 - 6
        min_k = minKmer - 2;
        max_k = maxKmer - 5;
    }else if(w < 1000){         // 5 - 9
        min_k = minKmer;
        max_k = maxKmer - 2;
    }else{                      // 5 - 11
        min_k = minKmer;
        max_k = maxKmer;   
    }
    
    int max_matches = -1;
    clear_rr(tmp_rr);  // clear the space for the result
    for(int k = min_k; k <= max_k; k++){
        clear_rr(rr);
        init_inputString(k, inputLen);
        find_tandem_repeat_sub(query_start, query_end, readID, inputLen, k, rr);
        
        if( max_matches < rr->Num_matches &&
           MIN_MATCH_RATIO < (float)rr->Num_matches/rr->repeat_len &&
           MIN_NUM_FREQ_UNIT < rr->Num_freq_unit &&
           MIN_PERIOD <= rr->rep_period )
        {
            max_matches = rr->Num_matches;
            *tmp_rr = *rr;
        }
    }
    *rr = *tmp_rr;
}

void insert_an_alignment(repeat_in_read rr){
    insert_an_alignment_into_set(
                                 rr.readID,
                                 rr.inputLen,
                                 rr.rep_start,
                                 rr.rep_end,
                                 rr.repeat_len,
                                 rr.rep_period,
                                 rr.Num_freq_unit,
                                 rr.Num_matches,
                                 rr.Num_mismatches,
                                 rr.Num_insertions,
                                 rr.Num_deletions,
                                 rr.Kmer,
                                 rr.ConsensusMethod,
                                 rr.string
                                 );
}

void remove_redundant_ranges_from_directional_index(int query_start, int query_end){
    for(int i=query_start; i<query_end; i++){
        if( directional_index[i] != -1 &&
           directional_index_end[i] < query_end )
        {
            directional_index[i]    = -1;
            directional_index_end[i]= -1;
            directional_index_w[i]  = -1;
        }
    }
}


void handle_one_TR(char *readID, int inputLen, int print_multiple_TR){

    struct timeval s_time, e_time;
    
    //
    // Locate overlapping regions of tandem repeats
    //
    gettimeofday(&s_time, NULL);
    
    int random_string_length;
    if(inputLen < MAX_WINDOW * 2){
        random_string_length = inputLen;
    }else{
        random_string_length = MAX_WINDOW * 2;
    }
    int DI_array_length = inputLen + random_string_length*2;
    
    fill_directional_index_with_end(DI_array_length, inputLen, random_string_length);

    gettimeofday(&e_time, NULL);
    time_range += (e_time.tv_sec - s_time.tv_sec) + (e_time.tv_usec - s_time.tv_usec)*1.0E-6;

#ifdef DEBUG_finding_ranges
    fprintf(stderr, "before removing non proper ranges --------------\n");
    for(int i=0; i<inputLen; i++){
        if(MIN_MAX_DI < directional_index[i] ){
            fprintf(stderr, "%i\t%i\t%i\t%i\t%f\n",
                    i,
                    directional_index_end[i],
                    (directional_index_end[i]- i),
                    directional_index_w[i],
                    directional_index[i]);
        }
    }
#endif
    
    //
    // Search for positions with tandem repeats
    //
    gettimeofday(&s_time, NULL);
    
    for(int query_start=0; query_start < inputLen; query_start++){
        int query_end = directional_index_end[query_start];
        int width         = directional_index_w[query_start];
        if( MIN_MAX_DI < directional_index[query_start] ){
            if(-1 < query_end){
                // Move onto de Bruijn graph construction
                find_tandem_repeat( query_start, query_end, width, readID, inputLen, &RRs[0], &RRs[1]);
                // Examine if a qualified TR is found
                if( RRs[0].repeat_len > 0 &&
                    RRs[0].rep_start + MIN_PERIOD * MIN_NUM_FREQ_UNIT < RRs[0].rep_end )
                {
                    insert_an_alignment(RRs[0]);
                    remove_redundant_ranges_from_directional_index(RRs[0].rep_start, RRs[0].rep_end);
                }
            }
        }
    }
    if(print_multiple_TR)
        chaining();
    else
        search_max();

    gettimeofday(&e_time, NULL);
    time_period += (e_time.tv_sec - s_time.tv_sec) + (e_time.tv_usec - s_time.tv_usec)*1.0E-6;
    
#ifdef DEBUG_finding_ranges
    fprintf(stderr, "after removing non proper ranges --------------\n");
    for(int i=0; i<inputLen; i++){
        if(MIN_MAX_DI < directional_index[i] ){
            fprintf(stderr, "%i\t%i\t%i\t%i\t%f\n",
                    i,
                    directional_index_end[i],
                    (directional_index_end[i]- i),
                    directional_index_w[i],
                    directional_index[i]);
        }
    }
#endif
    
}

void handle_one_read(char *readID, int inputLen, int read_cnt, int print_multiple_TR)
{
    handle_one_TR(readID, inputLen, print_multiple_TR);
}


