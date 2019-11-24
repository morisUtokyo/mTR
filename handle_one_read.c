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
            default: fprintf(stderr, "handle_one_read: fatal input char %i\n", val%4); exit(EXIT_FAILURE);
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
        default: fprintf(stderr, "handle_one_read: fatal input char %i\n", val); exit(EXIT_FAILURE);
    }
    return(return_char);
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
    
    strcpy( rr->readID, readID);
    rr->inputLen = inputLen;
    rr->Kmer     = k;
    
    struct timeval s, e;
    gettimeofday(&s, NULL);
    
    search_De_Bruijn_graph(query_start, query_end, rr);
    
    gettimeofday(&e, NULL);
    time_search_De_Bruijn_graph += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
    if(rr->rep_period < MIN_PERIOD){
        return;
    }else if(rr->rep_period * (query_end - query_start + 1) > WrapDPsize){
        fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
        clear_rr(rr);
    }else{
        wrap_around_DP(query_start, query_end, rr);
        // Polish the repeat unit if it is short and its coverage is small.
        // Otherwise, the coverage would be large enough.
        int coverage = rr->repeat_len / rr->rep_period;
        if( 5 <= coverage && coverage <= 20 && 10 < rr->rep_period){
            gettimeofday(&s, NULL);
            
            polish_repeat(rr);
            revise_representative_unit(rr);
            //wrap_around_DP(query_start, query_end, rr);
            
            gettimeofday(&e, NULL);
            time_polish += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
        }
        
    }
}

void find_tandem_repeat(int query_start, int query_end, int w, char *readID, int inputLen, repeat_in_read *rr, repeat_in_read *tmp_rr ){
    
    // A heuristic rule for setting the k-mer range for a de Bruijn graph to an nearly optimal one
    int min_k, max_k;
    if(w < 100){
        min_k = minKmer - 3;    // 2
        max_k = maxKmer - 5;    // 6  The accuracy gets worse if this is set to 5.
    }else if(w < 1000){
        min_k = minKmer - 3;    // 2
        max_k = maxKmer - 2;    // 9
    }else{
        min_k = minKmer;        // 5
        max_k = maxKmer;        // 11
    }
    
    int max_matches = -1;
    clear_rr(tmp_rr);  // clear the space for the result
    for(int k = min_k; k <= max_k; k++){
        clear_rr(rr);
        find_tandem_repeat_sub(query_start, query_end, readID, inputLen, k, rr);
        
        if( max_matches < rr->Num_matches &&
           MIN_MATCH_RATIO < (float)rr->Num_matches/rr->repeat_len &&
           MIN_NUM_FREQ_UNIT < rr->Num_freq_unit &&
           MIN_PERIOD <= rr->rep_period )
        {
            max_matches = rr->Num_matches;
            assign_rr(tmp_rr, rr);
        }
    }
    assign_rr(rr, tmp_rr);
}

void insert_an_alignment(repeat_in_read rr){
    insert_an_alignment_into_set(
                                 rr.readID,
                                 rr.inputLen,
                                 rr.rep_start,
                                 rr.rep_end,
                                 rr.repeat_len,
                                 rr.rep_period,
                                 rr.predicted_rep_period,
                                 rr.Num_freq_unit,
                                 rr.Num_matches,
                                 rr.Num_mismatches,
                                 rr.Num_insertions,
                                 rr.Num_deletions,
                                 rr.Kmer,
                                 rr.ConsensusMethod,
                                 rr.string,
                                 rr.string_score
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

void handle_one_TR(char *readID, int inputLen, int print_multiple_TR, int print_alignment){
    //
    // Locate overlapping regions of tandem repeats
    //
    int random_string_length;
    if(inputLen < MAX_WINDOW * 2){
        random_string_length = inputLen;
    }else{
        random_string_length = MAX_WINDOW * 2;
    }
    int DI_array_length = inputLen + random_string_length*2;
    
    struct timeval s_time_range, e_time_range;
    gettimeofday(&s_time_range, NULL);
    
    fill_directional_index_with_end(DI_array_length, inputLen, random_string_length);
    
    gettimeofday(&e_time_range, NULL);
    time_range += (e_time_range.tv_sec - s_time_range.tv_sec) + (e_time_range.tv_usec - s_time_range.tv_usec)*1.0E-6;
    
    //
    // Search for positions with tandem repeats
    //
    struct timeval s_time, e_time;
    gettimeofday(&s_time, NULL);
    
    for(int query_start=0; query_start < inputLen; query_start++){
        int query_end = directional_index_end[query_start];
        if(-1 < query_end && query_end < inputLen && MIN_MAX_DI < directional_index[query_start])
        {
            // Move onto de Bruijn graph construction
            clear_rr(&RRs[0]); clear_rr(&RRs[1]);
            int width     = directional_index_w[query_start];
            
            find_tandem_repeat( query_start, query_end, width, readID, inputLen, &RRs[0], &RRs[1]);
            
            query_counter++;
            // Examine if a qualified TR is found
            if( RRs[0].repeat_len > 0 &&
               RRs[0].rep_start + MIN_PERIOD * MIN_NUM_FREQ_UNIT < RRs[0].rep_end )
            {
                insert_an_alignment(RRs[0]);
                remove_redundant_ranges_from_directional_index(RRs[0].rep_start, RRs[0].rep_end);
            }
        }
    }
    struct timeval s_time_chaining, e_time_chaining;
    gettimeofday(&s_time_chaining, NULL);
    
    if(print_multiple_TR)
        chaining(print_alignment);
    else
        search_max(print_alignment);
    
    gettimeofday(&e_time_chaining, NULL);
    time_chaining += (e_time_chaining.tv_sec - s_time_chaining.tv_sec) + (e_time_chaining.tv_usec - s_time_chaining.tv_usec)*1.0E-6;

    gettimeofday(&e_time, NULL);
    time_period += (e_time.tv_sec - s_time.tv_sec) + (e_time.tv_usec - s_time.tv_usec)*1.0E-6;
    
    
}

void handle_one_read(char *readID, int inputLen, int read_cnt, int print_multiple_TR, int print_alignment)
{
    handle_one_TR(readID, inputLen, print_multiple_TR, print_alignment);
}


