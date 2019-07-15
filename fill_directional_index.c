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
#include "MT.h"



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

void assign_rr(repeat_in_read *rr_a, repeat_in_read *rr_b){
    rr_a->ID                = rr_b->ID;
    strcpy( rr_a->readID, rr_b->readID);
    rr_a->inputLen          = rr_b->inputLen;
    rr_a->rep_start         = rr_b->rep_start;
    rr_a->rep_end           = rr_b->rep_end;
    rr_a->repeat_len        = rr_b->repeat_len;
    rr_a->rep_period        = rr_b->rep_period;
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
    
    for(int i=0; i<inputLen + random_string_length*4 && i<MAX_INPUT_LENGTH; i++){ inputString_w_rand[i] = (int)(genrand_int32()%4);
    }
    
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
    // Compute the initial values for s(sum), q(squared sum), ip(inner product), and d(Manhattan distance)
    double s_0, s_1, s_2, q_0, q_1, q_2, ip_01, ip_12, d_01, d_12;
    s_0 = 0; s_1 = 0; s_2 = 0; q_0 = 0; q_1 = 0; q_2 = 0; ip_01 = 0; ip_12 = 0; d_01 = 0; d_12 = 0;
    for(int i=0; i<pow4[k]; i++){
        s_0    += vector0[i];
        s_1    += vector1[i];
        s_2    += vector2[i];
        q_0    += vector0[i] * vector0[i];
        q_1    += vector1[i] * vector1[i];
        q_2    += vector2[i] * vector2[i];
        ip_01  += vector0[i] * vector1[i];
        ip_12  += vector1[i] * vector2[i];
        d_01   += DIFF( vector0[i], vector1[i] );
        d_12   += DIFF( vector1[i], vector2[i] );
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
        double DI;
#ifdef Manhattan_Distance
        DI = (1 - d_12/s_1) - (1 - d_01/s_0);
        // d_01 and d_12 are the Manhattan distances, s_0 and s_1 are the sums of vector elements,
        // d_12/s_1 and d_01/s_0 are the mismatch ratios between the two vectors, and
        // 1 - d_12/s_1 and 1 - d_01/s_0 are the match ratios.
#else
        DI =  P_12 - P_01;
#endif
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
        d_01  -= DIFF( vector0[ dec_i_v0 ], vector1[ dec_i_v0 ] );
        // vector0, increment
        if(inc_i_v0 != dec_i_v0){
            q_0   -= pow(vector0[ inc_i_v0 ], 2);
            ip_01 -= vector0[ inc_i_v0 ] * vector1[ inc_i_v0 ];
            d_01  -= DIFF( vector0[ inc_i_v0 ], vector1[ inc_i_v0 ] );
        }
        // vector1, decrement
        q_1   -= pow(vector1[ dec_i_v1 ], 2);
        ip_12 -= vector1[ dec_i_v1 ] * vector2[ dec_i_v1 ];
        d_12  -= DIFF( vector1[ dec_i_v1 ], vector2[ dec_i_v1 ] );
        // vector 1, increment
        if(inc_i_v1 != dec_i_v1){
            q_1   -= pow(vector1[ inc_i_v1 ], 2);
            if(inc_i_v1 != dec_i_v0){
                ip_01 -= vector0[ inc_i_v1 ] * vector1[ inc_i_v1 ];
                d_01  -= DIFF( vector0[ inc_i_v1 ], vector1[ inc_i_v1 ] );
            }
            ip_12 -= vector1[ inc_i_v1 ] * vector2[ inc_i_v1 ];
            d_12  -= DIFF( vector1[ inc_i_v1 ], vector2[ inc_i_v1 ] );
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
                d_12  -= DIFF( vector1[ inc_i_v2 ], vector2[ inc_i_v2 ] );
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
        d_01  += DIFF( vector0[ dec_i_v0 ], vector1[ dec_i_v0 ] );
        // vector0, increment
        if(inc_i_v0 != dec_i_v0){
            q_0   += pow(vector0[ inc_i_v0 ], 2);
            ip_01 += vector0[ inc_i_v0 ] * vector1[ inc_i_v0 ];
            d_01  += DIFF( vector0[ inc_i_v0 ], vector1[ inc_i_v0 ] );
        }
        // vector1, decrement
        q_1   += pow(vector1[ dec_i_v1 ], 2);
        // ip_01 += vector0[ dec_i_v1 ] * vector1[ dec_i_v1 ];
        // duplicated as inc_i_v0 = dec_i_v1
        ip_12 += vector1[ dec_i_v1 ] * vector2[ dec_i_v1 ];
        d_12  += DIFF( vector1[ dec_i_v1 ], vector2[ dec_i_v1 ] );
        // vector 1, increment
        if(inc_i_v1 != dec_i_v1){
            q_1   += pow(vector1[ inc_i_v1 ], 2);
            if(inc_i_v1 != dec_i_v0){
                ip_01 += vector0[ inc_i_v1 ] * vector1[ inc_i_v1 ];
                d_01  += DIFF( vector0[ inc_i_v1 ], vector1[ inc_i_v1 ] );
            }
            ip_12 += vector1[ inc_i_v1 ] * vector2[ inc_i_v1 ];
            d_12  += DIFF( vector1[ inc_i_v1 ], vector2[ inc_i_v1 ] );
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
                d_12  += DIFF( vector1[ inc_i_v2 ], vector2[ inc_i_v2 ] );
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
            // Search for a local minimum as the end position
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


void fill_directional_index_with_end(int DI_array_length, int inputLen, int random_string_length){
    
    struct timeval s_time, e_time;
    gettimeofday(&s_time, NULL);
    
    // initialize directional _index
    for(int i=0; i<DI_array_length; i++){
        directional_index[i]   = -1;
        directional_index_end[i] = -1;
    }
    int k;
    
    // For finding short tandem repeats of frequency <= 10 and of length 2-10
    k = 3;
    init_inputString_surrounded_by_random_seq(k, inputLen, random_string_length);
        // Put random sequences of the input length before and after the input string
    for(int w = 20; w < 50  && w < inputLen/2; w += 10) // w = 20, 30, 40
    {
        fill_directional_index_tmp( DI_array_length, w, k, inputLen, random_string_length);
        put_local_maximum_into_directional_index( DI_array_length, w );
    }

    // For finding longer tandem repeats
    k = 5;
    init_inputString_surrounded_by_random_seq(k, inputLen, random_string_length);
        // Put random sequences of the input length before and after the input string
    for(int w = MIN_WINDOW; w < MAX_WINDOW && w < inputLen/2; )
    {
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
    
    gettimeofday(&e_time, NULL);
    time_range += (e_time.tv_sec - s_time.tv_sec) + (e_time.tv_usec - s_time.tv_usec)*1.0E-6;
}
