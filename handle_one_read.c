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
#include "MT.h"
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


void clear_rr(repeat_in_read *rr_a){
    rr_a->ID                = -1;
    strcpy( rr_a->readID, "");
    rr_a->inputLen          = -1;
    rr_a->rep_start         = -1;
    rr_a->rep_end           = -1;
    rr_a->repeat_len        = -1;
    rr_a->rep_period        = -1;
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

void find_tandem_repeat(int query_start, int query_end, int w, char *readID, int inputLen, int k, repeat_in_read *rr, repeat_in_read *tmp_rr ){
    
    // A heuristic rule for setting the maximum k-mer to an oprimal value
    int min_k = minKmer;

    int diff;
    if(w < 100){
        diff = 5;
    }else if(w < 1000){
        diff = 2;
    }else{
        diff=0;
    }
    int max_k = maxKmer - diff;
    
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

double DI_index(int *vector0, int *vector1, int *vector2, int k){
    double s_0 = 0;
    double s_1 = 0;
    double s_2 = 0;
    double q_0 = 0;
    double q_1 = 0;
    double q_2 = 0;
    double ip_01 = 0;
    double ip_12 = 0;
    for(int i=0; i<pow4[k]; i++){
        s_0    += vector0[i];
        s_1    += vector1[i];
        s_2    += vector2[i];
        q_0    += vector0[i] * vector0[i];
        q_1    += vector1[i] * vector1[i];
        q_2    += vector2[i] * vector2[i];
        ip_01 += vector0[i] * vector1[i];
        ip_12 += vector1[i] * vector2[i];
    }

#ifdef DEBUG_incremental
    fprintf(stderr, "Batch\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", s_0, s_1, s_2, q_0, q_1, q_2, ip_01, ip_12);
#endif
    
    double sd_0    = sqrt(q_0 * pow4[k] - s_0 * s_0);
    double sd_1    = sqrt(q_1 * pow4[k] - s_1 * s_1);
    double sd_2    = sqrt(q_2 * pow4[k] - s_2 * s_2);
    
    double P_01, P_12;
    if(sd_0*sd_1 > 0){
        P_01 = (ip_01 * pow4[k] - s_0 * s_1)/(sd_0 * sd_1);
    }else{
        P_01 = 0;
    }
    if(sd_1*sd_2 > 0){
        P_12 = (ip_12 * pow4[k] - s_1 * s_2)/(sd_1 * sd_2);
    }else{
        P_12 = 0;
    }
    double DI = P_12 - P_01;
    return(DI);
}

void init_inputString_surrounded_by_random_seq(int k, int inputLen, int random_string_length){
    init_genrand(0);
    
    for(int i=0; i<inputLen + random_string_length*4 && i<MAX_INPUT_LENGTH; i++){ inputString_w_rand[i] = (int)(genrand_int32()%4);}
    
    for(int i = 0; i < random_string_length; i++){
        inputString_w_rand[i] = (int)(genrand_int32()%4);
    }
    for(int i = 0; i < inputLen; i++){
        inputString_w_rand[i + random_string_length] = orgInputString[i];
    }
    for(int i = 0; i < random_string_length; i++){
        inputString_w_rand[i + random_string_length + inputLen] = (int)(genrand_int32()%4);
    }
    int tmp = 0;
    for(int i=0; i<(k-1); i++){
        tmp = 4 * tmp + inputString_w_rand[i];  // compute 4 decimal of first k-1 letters
    }
    for(int i=0; i<(inputLen + random_string_length*2 - k+1); i++){
        inputString_w_rand[i] = 4 * tmp + inputString_w_rand[i+k-1];
        tmp = inputString_w_rand[i] % pow4[k-1]; //　remainder, compute 4 decimal of length k-1
        if(tmp < 0){
            fprintf(stderr, "fatal error at %i\t %i \t%i\n", i, inputString_w_rand[i], pow4[k-1] );
            exit(EXIT_FAILURE);
        }
    }
}

void fill_directional_index_tmp(int DI_array_length, int w, int k, int inputLen, int random_string_length){
    // We use inputLen and random_string_length for analyzing patterns of DI and Pearson's CC only.
    
    for(int i=0; i<DI_array_length; i++){
        directional_index_tmp[i] = -1;
    }
    
    // initialize vectors
    for(int i=0; i < 4 * BLK; i++){
        vector0[i] = 0;
        vector1[i] = 0;
        vector2[i] = 0;
    }

    for(int i=0; i<w; i++){
        vector0[ inputString_w_rand[i] ]++;
        vector1[ inputString_w_rand[i + w] ]++;
        vector2[ inputString_w_rand[i + w*2] ]++;
    }
    // https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    // Compute the initial values for s(sum), q(squared sum), and ip(inner product)
    double s_0, s_1, s_2, q_0, q_1, q_2, ip_01, ip_12;
    s_0 = 0; s_1 = 0; s_2 = 0; q_0 = 0; q_1 = 0; q_2 = 0; ip_01 = 0; ip_12 = 0;
    for(int i=0; i<pow4[k]; i++){
        s_0    += vector0[i];
        s_1    += vector1[i];
        s_2    += vector2[i];
        q_0    += vector0[i] * vector0[i];
        q_1    += vector1[i] * vector1[i];
        q_2    += vector2[i] * vector2[i];
        ip_01  += vector0[i] * vector1[i];
        ip_12  += vector1[i] * vector2[i];
    }
    // Shift vector0/1/2 starting from i=1 can calculate directional indexes
    for(int i=0; i < (DI_array_length - w - random_string_length - k + 1); i++){
        if(s_0 == s_1 && s_1 == s_2){}
        else{
            fprintf(stderr, "inconsistency in fill_directional_index_tmp\n");
            exit(EXIT_FAILURE);
        }
        // Compute the directional index at i
        double sd_0    = sqrt(q_0 * pow4[k] - s_0 * s_0);
        double sd_1    = sqrt(q_1 * pow4[k] - s_1 * s_1);
        double sd_2    = sqrt(q_2 * pow4[k] - s_2 * s_2);
        double P_01, P_12;
        if(sd_0*sd_1 > 0){
            P_01 = (ip_01 * pow4[k] - s_0 * s_1)/(sd_0 * sd_1);
        }else{  // When one of the two SDs is zero.
            P_01 = 0;
        }
        if(sd_1*sd_2 > 0){
            P_12 = (ip_12 * pow4[k] - s_1 * s_2)/(sd_1 * sd_2);
        }else{
            P_12 = 0;
        }
        double DI =  P_12 - P_01;
        // Note that DI is computed for position i+w (NOT i !!)
        directional_index_tmp[i+w] = DI;
        
#ifdef DUMP_DI_PCC
        int real_pos = i+w-random_string_length;
        if( 0 <= real_pos && real_pos < inputLen){
            fprintf(stdout, "%i,%f\n", real_pos, P_12);     // Pearson's correlation coefficients
            //fprintf(stdout, "%i,%f\n", real_pos, DI);   // Directional index
        }
#endif
        // Incremental updates of s, q, and ip
        // Determine positions to be updated
        // dec_i_v0 means the position in vector0 was decremented.
        int dec_i_v0, inc_i_v0, dec_i_v1, inc_i_v1, dec_i_v2, inc_i_v2;
        dec_i_v0 = inputString_w_rand[i];
        inc_i_v0 = inputString_w_rand[i + w];
        dec_i_v1 = inputString_w_rand[i + w];
        inc_i_v1 = inputString_w_rand[i + w*2];
        dec_i_v2 = inputString_w_rand[i + w*2];
        inc_i_v2 = inputString_w_rand[i + w*3];
        
        // Subtract changes before updates
        // vector0, decrement
        q_0   -= pow(vector0[ dec_i_v0 ], 2);
        ip_01 -= vector0[ dec_i_v0 ] * vector1[ dec_i_v0 ];
        // vector0, increment
        if(inc_i_v0 != dec_i_v0){
            q_0   -= pow(vector0[ inc_i_v0 ], 2);
            ip_01 -= vector0[ inc_i_v0 ] * vector1[ inc_i_v0 ];
        }
        // vector1, decrement
        q_1   -= pow(vector1[ dec_i_v1 ], 2);
            // ip_01 -= vector0[ dec_i_v1 ] * vector1[ dec_i_v1 ];
            // duplicated as inc_i_v0 = dec_i_v1
        ip_12 -= vector1[ dec_i_v1 ] * vector2[ dec_i_v1 ];
        // vector 1, increment
        if(inc_i_v1 != dec_i_v1){
            q_1   -= pow(vector1[ inc_i_v1 ], 2);
            if(inc_i_v1 != dec_i_v0){
                ip_01 -= vector0[ inc_i_v1 ] * vector1[ inc_i_v1 ];
            }
            ip_12 -= vector1[ inc_i_v1 ] * vector2[ inc_i_v1 ];
        }
        // vector2, decrement
        q_2   -= pow(vector2[ dec_i_v2 ], 2);
            // ip_12 -= vector1[ dec_i_v2 ] * vector2[ dec_i_v2 ];
            // duplicated as inc_i_c1 == dec_i_v2
        // vector2, increment
        if(inc_i_v2 != dec_i_v2){
            q_2   -= pow(vector2[ inc_i_v2 ], 2);
            if(inc_i_v2 != dec_i_v1){
                ip_12 -= vector1[ inc_i_v2 ] * vector2[ inc_i_v2 ];
            }
        }
        
        // Update vectors and save the previous values before decrement/increment vectors
        vector0[ dec_i_v0 ]--;  s_0--;
        vector0[ inc_i_v0 ]++;  s_0++;
        vector1[ dec_i_v1 ]--;  s_1--;
        vector1[ inc_i_v1 ]++;  s_1++;
        vector2[ dec_i_v2 ]--;  s_2--;
        vector2[ inc_i_v2 ]++;  s_2++;
        // Add changes after updates
        // vector0, decrement
        q_0   += pow(vector0[ dec_i_v0 ], 2);
        ip_01 += vector0[ dec_i_v0 ] * vector1[ dec_i_v0 ];
        // vector0, increment
        if(inc_i_v0 != dec_i_v0){
            q_0   += pow(vector0[ inc_i_v0 ], 2);
            ip_01 += vector0[ inc_i_v0 ] * vector1[ inc_i_v0 ];
        }
        // vector1, decrement
        q_1   += pow(vector1[ dec_i_v1 ], 2);
            // ip_01 += vector0[ dec_i_v1 ] * vector1[ dec_i_v1 ];
            // duplicated as inc_i_v0 = dec_i_v1
        ip_12 += vector1[ dec_i_v1 ] * vector2[ dec_i_v1 ];
        // vector 1, increment
        if(inc_i_v1 != dec_i_v1){
            q_1   += pow(vector1[ inc_i_v1 ], 2);
            if(inc_i_v1 != dec_i_v0){
                ip_01 += vector0[ inc_i_v1 ] * vector1[ inc_i_v1 ];
            }
            ip_12 += vector1[ inc_i_v1 ] * vector2[ inc_i_v1 ];
        }
        // vector2, decrement
        q_2   += pow(vector2[ dec_i_v2 ], 2);
            // ip_12 += vector1[ dec_i_v2 ] * vector2[ dec_i_v2 ];
            // duplicated as inc_i_c1 == dec_i_v2
        // vector2, increment
        if(inc_i_v2 != dec_i_v2){
            q_2   += pow(vector2[ inc_i_v2 ], 2);
            if(inc_i_v2 != dec_i_v1){
                ip_12 += vector1[ inc_i_v2 ] * vector2[ inc_i_v2 ];
            }
        }
    }
}

void put_local_maximum_into_directional_index(int DI_array_length, int w){
    // Search for local maximums
    double local_max = -1;  // Set it to the mimimum -1
    int   local_max_i = -1;
    for(int i=0; i < DI_array_length; i++){
        if( local_max < directional_index_tmp[i] ){
            local_max = directional_index_tmp[i];
            local_max_i = i;
        }
        if(local_max_i + w < i &&
           directional_index[local_max_i] < local_max &&
           MIN_MAX_DI < local_max )
        {
            // The position, local_max_i, was updated more than w before the current position i.
            // It must be also greater than the maximum at the position.
            // Search for a local minimum
            double local_min = 1;   // Set it to the maximum 1
            int local_min_j = local_max_i;
            for(int j=local_max_i; j < DI_array_length; j++){
                // Search the temporary directional index array !!
                if( local_min > directional_index_tmp[j] ){
                    local_min = directional_index_tmp[j];
                    local_min_j = j;
                }
                if(local_min_j + w < j){
                    directional_index[local_max_i] = local_max;
                    directional_index_w[local_max_i] = w;
                    directional_index_end[local_max_i] = local_min_j + w;  // set the end position
                    i = local_min_j + w;
                    break;
                }
            }
            local_max = -1;
            local_min = 1;
        }
    }
}


void fill_directional_index_with_end(int DI_array_length, int inputLen, int k, int random_string_length){
    
    // Put random sequences of the input length before and after the input string
    init_inputString_surrounded_by_random_seq(k, inputLen, random_string_length);
    
    // initialize directional _index
    for(int i=0; i<DI_array_length; i++){
        directional_index[i]   = -1;
        directional_index_end[i] = -1;
    }
    for(int w = min_window_size; w < max_window_size && w < inputLen/2; ){
        fill_directional_index_tmp( DI_array_length, w, k, inputLen, random_string_length);
        put_local_maximum_into_directional_index( DI_array_length, w );
        // Sizes of windows
        if(w < 100){
            w += 20;
        }else if(w < 1000){
            w += 300;
        }else if(w < 20000){
            w += 2000;
        }else{
            w += 10000;
        }
    }
    
    //  By removing random sequences of both ends, retain directional indexes of the input
    for(int i=0; i<inputLen; i++){
        int shifted_i = i + random_string_length;
        directional_index[i]    = directional_index[shifted_i];
        directional_index_end[i]= directional_index_end[shifted_i] - random_string_length;
        directional_index_w[i]  = directional_index_w[shifted_i];
    }
    for(int i=inputLen; i<DI_array_length; i++){
        directional_index[i]    = -1;
        directional_index_end[i]= -1;
        directional_index_w[i]  = -1;
    }
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
    
    int k = 4;  // Setting k to 3 is inferior to k = 4 when units are of length  5.
    int random_string_length;
    if(inputLen < max_window_size * 2){
        random_string_length = inputLen;
    }else{
        random_string_length = max_window_size * 2;
    }
    int DI_array_length = inputLen + random_string_length*2;
    fill_directional_index_with_end(DI_array_length, inputLen, k, random_string_length);

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
                find_tandem_repeat( query_start, query_end, width, readID, inputLen, k, &RRs[0], &RRs[1]);
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


