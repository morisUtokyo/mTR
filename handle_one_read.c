//
//  handle_one_read.c
//  
//
//  Created by Shinichi Morishita
//

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

void find_tandem_repeat_sub(int max_start, int max_end, char *readID, int inputLen, int Kmer, repeat_in_read *rr ){
    
    //---------------------------------------------------------------------------
    // In the range [max_start, max_end],
    // compute the representative unit string by using a progressive multiple alignment or
    // by traversing the De Bruijn graph of all k-mers in a greedy manner
    //---------------------------------------------------------------------------
    
    // Sort k-mers using counting sort and store the positions of k-mers in sortedString
    int pow4k = 1;    // 4^{k}  e.g., 4^(4) = 256
    for(int i=0; i<Kmer; i++){ pow4k = 4 * pow4k;}
    int pow4k_1 = pow4k/4;
    for(int i = 0; i < pow4k; i++){ count[i] = 0;}  // Initialization
    for(int i = max_start; i <= max_end; i++){ count[ inputString[i] ]++; }    // Perform counting
    for(int i = 1; i < pow4k; i++){ count[i] = count[i-1] + count[i]; }
    for(int i=0; i<MAX_INPUT_LENGTH; i++){ sortedString[i] = 0; } // Initialization
    for(int i=max_end; max_start <= i; i--){
        sortedString[ --count[inputString[i]] ] = i;
    }
    
    // Initialization
    for(int i = 0; i < MAX_PERIOD; i++){
        count_period_all[i] = 0;
    }
    for(int j=0; j<MAX_INPUT_LENGTH; j++){
        freq_interval_len[j] = 0;
    }
    for(int i=max_start; i <= max_end; i++){
        int local_start = count[inputString[i]];
        int local_end;
        if(inputString[i] == pow4k-1){
            local_end = inputLen;
        }else{
            local_end = count[ inputString[i]+1 ];
        }
        freq_interval_len[i] = (double)(local_end - local_start);
        
        for(int j = local_start; j < (local_end - 1); j++){
            int aPeriod = abs(sortedString[j+1] - sortedString[j]);
            if(aPeriod < MAX_PERIOD){
                count_period_all[ aPeriod ]++;
            }
        }
    }
    //　Locate the initial position (k-mer) with the maximum frequency
    int max_pos = max_start;
    int max_freq   = (int)freq_interval_len[max_pos];
    for(int j=max_start; j <= max_end; j++){
        if(max_freq < freq_interval_len[j]){
            max_pos  = j;
            max_freq = (int)freq_interval_len[max_pos];
        }
    }
    // Compute repeat period with the maximum frequency
    int rep_period = 2;
    max_freq = 0;
    for(int p = 2; p < MAX_PERIOD; p++){
        if(max_freq < count_period_all[p]){
            rep_period = p;
            max_freq = count_period_all[p];
        }
    }
    // Initialization
    int actual_rep_period;
    for(int i = 0; i < MAX_PERIOD; i++){
        rep_unit_string[i] = 0;
    }
    
    // Traverse the De Bruijn graph of all k-mers in a greedy manner
    actual_rep_period =
    search_De_Bruijn_graph(max_start, max_end, max_pos, rep_period, inputLen, pow4k_1);
    
    int ConsensusMethod;
    //  If a repeat unit is found, actual_rep_period > 0, and = 0 otherwise.
    if(actual_rep_period > 0){
        // A De Bruijin graph search is successful.
        ConsensusMethod = DeBruijnGraphSearch;
    }else{
        // If the De Bruijn graph search fails, try a progressive  multiple alignment.
        ConsensusMethod = ProgressiveMultipleAlignment;
        actual_rep_period = progressive_multiple_alignment( max_start, max_end, max_pos, rep_period, Kmer, inputLen, pow4k);
    }
    
    //---------------------------------------------------------------------------
    // Compute the accuracy of the representative unit string by wrap-around DP
    //---------------------------------------------------------------------------
    
    if(actual_rep_period * (max_end - max_start + 1) > WrapDPsize){
        fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
        clear_rr(rr);
    }else if(actual_rep_period > 0){
        int actual_repeat_len, Num_freq_unit, Num_matches, Num_mismatches, Num_insertions, Num_deletions;
        
        wrap_around_DP(rep_unit_string,
                       actual_rep_period,
                       &orgInputString[max_start],
                       (max_end - max_start + 1),
                       &actual_repeat_len,  &Num_freq_unit,
                       &Num_matches,        &Num_mismatches,
                       &Num_insertions,     &Num_deletions);
        
        strcpy( rr->readID, readID);
        rr->inputLen           = inputLen;
        rr->max_start          = max_start;
        rr->max_end            = max_end;
        rr->rep_period         = rep_period;
        rr->actual_rep_period  = actual_rep_period;
        rr->Kmer      = Kmer;
        rr->Num_freq_unit      = 0;
        rr->ConsensusMethod    = ConsensusMethod;
        rr->Num_freq_unit      = Num_freq_unit;
        rr->Num_matches        = Num_matches;
        rr->Num_mismatches     = Num_mismatches;
        rr->Num_insertions     = Num_insertions;
        rr->Num_deletions      = Num_deletions;
        
        rr->actual_repeat_len = actual_repeat_len;
        print_4_decimal_array(rep_unit_string, actual_rep_period, rr->string);
        freq_2mer_array(rep_unit_string, actual_rep_period, rr->freq_2mer);
    }
}

void init_inputString(int k, int inputLen){
    // Encode the raw input string into 4 decimals of length k
    int pow4k, pow4k_1;
    pow4k_1 = 1;    // 4^{k-1}  e.g., 4^(4-1) = 64
    for(int i=0; i<(k-1); i++){
        pow4k_1 = 4 * pow4k_1;
    }
    for(int i=0; i<inputLen; i++){
        inputString[i] = orgInputString[i];
    }
    int tmp = 0;
    for(int i=0; i<(k-1); i++){
        tmp = 4 * tmp + inputString[i];  // compute 4 decimal of first k-1 letters
    }
    for(int i=0; i<(inputLen-k+1); i++){
        inputString[i] = 4 * tmp + inputString[i+k-1];
        tmp = inputString[i] % pow4k_1; //　remainder, compute 4 decimal of length k-1
        if(tmp < 0){
            fprintf(stderr, "fatal error at %i\t %i \t%i\n", i, inputString[i], pow4k_1);
            exit(EXIT_FAILURE);
        }
    }
}


void find_tandem_repeat(int max_start, int max_end, char *readID, int inputLen,  repeat_in_read *rr, repeat_in_read *tmp_rr ){
    
    int max_matches = -1;
    clear_rr(tmp_rr);  // clear the space for the result
    
    for(int k = minKmer; k <= maxKmer; k++){
        clear_rr(rr);
        init_inputString(k, inputLen);
        find_tandem_repeat_sub(max_start, max_end, readID, inputLen, k, rr);
        
        if( max_matches < rr->Num_matches &&
            MIN_MATCH_RATIO < (float)rr->Num_matches/rr->actual_repeat_len &&
            MIN_NUM_FREQ_UNIT < rr->Num_freq_unit )
        {
            max_matches = rr->Num_matches;
            *tmp_rr = *rr;
        }
    }
    *rr = *tmp_rr;
}


double DI_index(int *v0, int *v1, int *v2, int k){
    int pow4k = 1;
    for(int i=0; i<k; i++){
        pow4k = 4 * pow4k;
    }
    double s_0 = 0;
    double s_1 = 0;
    double s_2 = 0;
    double q_0 = 0;
    double q_1 = 0;
    double q_2 = 0;
    double cov_01 = 0;
    double cov_12 = 0;
    for(int i=0; i<pow4k; i++){
        s_0    += v0[i];
        s_1    += v1[i];
        s_2    += v2[i];
        q_0    += v0[i] * v0[i];
        q_1    += v1[i] * v1[i];
        q_2    += v2[i] * v2[i];
        cov_01 += v0[i] * v1[i];
        cov_12 += v1[i] * v2[i];
    }
    double sd_0    = sqrt(q_0 * pow4k - s_0 * s_0);
    double sd_1    = sqrt(q_1 * pow4k - s_1 * s_1);
    double sd_2    = sqrt(q_2 * pow4k - s_2 * s_2);
    
    // https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#For_a_sample
    double P_01, P_12;
    if(sd_0*sd_1 > 0){
        P_01 = (cov_01 * pow4k - s_0 * s_1)/(sd_0 * sd_1);
    }else{
        P_01 = 0;
    }
    if(sd_1*sd_2 > 0){
        P_12 = (cov_12 * pow4k - s_1 * s_2)/(sd_1 * sd_2);
    }else{
        P_12 = 0;
    }
    double DI = P_12 - P_01;
    return(DI);
}

void fill_directional_index(int inputLen, int w, int k){
    // initialized inputString with a given value for k
    init_inputString(k, inputLen);
    // initialize directional _index
    for(int i=0; i<MAX_INPUT_LENGTH; i++){
        directional_index[i] = 0;
    }
    // initialize vectors
    for(int i=0; i < 4 * BLK; i++){
        vector0[i] = 0; vector1[i] = 0; vector2[i] = 0;
    }
    // Assume that vector0, the window ending at (i-1), is the empty vector filled with 0
    for(int i=0; i < w && i < (inputLen - w - k + 1); i++){
        vector1[ inputString[i] ]++;
        vector2[ inputString[i + w] ]++;
    }
    // Shift vector0/1/2 starting from i=1 can calculate directional indexes
    for(int i=0; i < (inputLen - w - k + 1); i++){
        if(0 <= i-w){
            vector0[ inputString[i-w] ]--;
        }
        vector0[ inputString[i] ]++;
        vector1[ inputString[i] ]--;
        vector1[ inputString[i + w] ]++;
        vector2[ inputString[i + w] ]--;
        if(w*2 < inputLen - k + 1){
            vector2[ inputString[i + w*2] ]++;
        }
        // starting index at the i-th position
        directional_index[i] = (double)DI_index(vector0, vector1, vector2, k);
        //printf("%c\t%0.2f ", dec2char(orgInputString[i]), directional_index[i]);
    }
}

int find_best_candidate_region(int inputLen, int w, int k, int search_pos, int *max_start, int *max_end, double *max_DI_answer, int print_multiple_TR){
    
    int found = 0;
    int max_pos = search_pos;
    double max_DI = 0;
    int min_pos = search_pos;
    double min_DI = 0;
    for(int i=search_pos; i < (inputLen - w - k + 1); i++){
        if(max_DI < directional_index[i]){
            max_pos = i;
            max_DI  = directional_index[max_pos];
        }
        if(directional_index[i] < min_DI){
            min_pos = i;
            min_DI  = directional_index[min_pos];
        }
        if( print_multiple_TR == 1 &&
            max_pos + w < i &&
            max_pos + w < min_pos &&
            min_pos + w < i )
        {   // max_pos and min_pos are locally maximum and minimum, respectively and are distant apart >= w bases.
            found = 1;
            break;
        }
    }
    if(print_multiple_TR == 0 && max_pos < min_pos){
        // For the mode of finding one tandem repeat print_multiple_TR == 0)
        found = 1;
    }
    if(found == 1){
        *max_start = max_pos;
        *max_end   = min_pos + w;
        *max_DI_answer = max_DI;
    }
    return(found);
}

void handle_one_read(char *readID, int inputLen, int read_cnt, int print_multiple_TR){
    
    int search_pos = 0;
    while(search_pos < inputLen){
        double max_measure = 0;
        int found = 0;
        int max_w, max_k, tmp_start, tmp_end, max_start, max_end;
        double tmp_DI = 0;
        double max_DI = 0;
        
        int w0=64;
        for(int k=3; k<=5; k++){ // w0 = 4^k
            int w, iter;
            if(k == 3){         //  w = 32, 64 or 128.  w/4^k = 1/2, 1 or 2
                w = w0/2;   iter = 3;
            }else if(k ==4){    //  w = 256, 512.  w/4^k = 1 or 2
                w = w0;     iter = 2;
            }else{              // w = 1024, 2048, w/4^k = 1 or 2
                w = w0;     iter = 2;
            }
            for(int i=0; i<iter; i++){
                fill_directional_index(inputLen, w, k);
                int found_one = find_best_candidate_region(inputLen, w, k, search_pos, &tmp_start, &tmp_end, &tmp_DI, print_multiple_TR);
                
                double tmp_measure;
                tmp_measure = tmp_DI * (1 + 0.2 * (k-3));
                //tmp_measure = log(max_end - max_start + 1);
                
                if(found_one == 1 && max_measure < tmp_measure) {
                    max_measure = tmp_measure;
                    max_w = w;
                    max_k = k;
                    max_start = tmp_start;
                    max_end   = tmp_end;
                    max_DI    = tmp_DI;
                    found = 1;
                }
                w = w * 2;
            }
            w0 = w0 * 4;
        }

        if(found){
#ifdef DEBUG_window_kmer
            printf("w=%i\tk=%i\t[%i, %i]\n", max_w, max_k, max_start, max_end);
#endif
            fill_directional_index(inputLen, max_w, max_k);
            find_tandem_repeat(max_start, max_end, readID, inputLen, &RRs[0], &RRs[1]);
            
            if(RRs[0].actual_repeat_len > 0){
                print_one_repeat_in_read(RRs[0]);
            }
            if(print_multiple_TR == 1){
                search_pos = max_end;
            }else{
                break;
            }
        }else{
            break;
        }

    }
}


