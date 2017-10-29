//
//  Consensus.c
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

int maxFreqKmer(int inputLen, int pow4k){
    int maxCnt = 0;
    int AnsKmer = 0;
    for(int i = 1; i < pow4k; i++){
        int Cnt;
        if(i == pow4k - 1){
            Cnt = inputLen - count[i];
        }else{
            Cnt = count[i+1] - count[i];
        }
        if(maxCnt < Cnt){
            maxCnt = Cnt;
            AnsKmer = i;
        }
    }
    return(AnsKmer);
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

void search_De_Bruijn_graph(int pow4k, repeat_in_read *tmp_rr, repeat_in_read *rr){
    
    int max_start   = tmp_rr->max_start;
    int max_end     = tmp_rr->max_end;
    int rep_period  = tmp_rr->rep_period;
    int Kmer        = tmp_rr->Kmer;
    int inputLen    = tmp_rr->inputLen;
    int actual_rep_period = 0;
    // Starting from the initial k-mer, traverse the De Bruijn graph of all k-mers in a greedy manner
    
    int Node = maxFreqKmer(inputLen, pow4k);
    
    int pow4k_1 = pow4k / 4;
    int start, end, tmp_count;
    int initialNode = Node;
    int nextNode;
    
    int rep_unit_string[MAX_PERIOD];
    for(int k=0; k< MAX_PERIOD; k++){
        rep_unit_string[k] = Node / pow4k_1;
        
        int max_lsd = 0;
        int max_count_lsd = 0;
        
        for(int lsd=0; lsd < 4; lsd++){ // lsd: least significant (4) digit
            nextNode = 4 *(Node % pow4k_1) + lsd;
            start = count[nextNode];
            if(nextNode == 4*pow4k_1 - 1){
                end = inputLen;
            }else{
                end = count[ nextNode+1 ];
            }
            
            tmp_count = end - start + 1;
            
            if(max_count_lsd < tmp_count){
                max_count_lsd = tmp_count;
                max_lsd = lsd;
            }
        }
        Node = 4 * (Node % pow4k_1) + max_lsd;
        
        // If you hit the initial node, terminate the search.
        if(Node == initialNode && k+1 < MAX_PERIOD){
            actual_rep_period = k+1;
            break;
        }
    }
    
    rr->actual_rep_period = actual_rep_period;
    rr->ConsensusMethod   = DeBruijnGraphSearch;
    if(actual_rep_period > 0){
        wrap_around_DP(rep_unit_string,
                       actual_rep_period,
                       &orgInputString[max_start],
                       (max_end - max_start + 1),
                       rr);
        print_4_decimal_array(rep_unit_string, actual_rep_period, rr->string);
        freq_2mer_array(rep_unit_string, actual_rep_period, rr->freq_2mer);
    }
}

void progressive_multiple_alignment(
           int pow4k, repeat_in_read *tmp_rr, repeat_in_read *rr){
    
    int max_start   = tmp_rr->max_start;
    int max_end     = tmp_rr->max_end;
    int rep_period  = tmp_rr->rep_period;
    int Kmer        = tmp_rr->Kmer;
    int inputLen    = tmp_rr->inputLen;

    // Start from the positions associated with k-mers of the maximum frequency in "count" data structure ranging from "start" to "end defined below
    
    int MaxFreqKmer = maxFreqKmer(inputLen, pow4k);
    
    int start, end;
    start = count[MaxFreqKmer];
    if(MaxFreqKmer == pow4k - 1){
        end = inputLen;
    }else{
        end = count[MaxFreqKmer + 1];
    }
    
#ifdef DEBUG_progressive_multiple_alignment
    printf("max_start = %i\tmax_end = %i\trep_period =%i\tstart=%i\tend=%i\tKmer = %i\n", max_start, max_end, rep_period, start, end, Kmer);
#endif
    
    
    int *matDP = WrapDP;  // Reuse WrapDP by renaming WrapDP by oneDP
    int next = rep_period + 1;
    
    int consensus[MAX_PERIOD][5];
    for(int i=0; i<MAX_PERIOD; i++){
        for(int j=0; j<5; j++){
            consensus[i][j] = 0;
        }
    }
    int gaps[MAX_PERIOD+1][4];
    
    int k = start;
    // Move k into the repeat region
    for( ; ( sortedString[k] < max_start) && (k < end); k++);
    int max_pos = sortedString[k];
    
    // A progressive multiple alignment between
    // the reference starting from max_pos and
    // one string starting from sortedString[k] in oneInputString
    
    for(; ( sortedString[k] <= max_end) && (k < end); k++){
        // Perform a global alignment
        // i scans the reference starting from max_pos in oneInputString.
        // j scans one string starting from sortedString[k] in oneInputString.
        
        int i, j;
        // Initialize the boundaries
        for(i = 0; i <= rep_period; i++){
            matDP[ i * next + 0] = INDEL_PENALTY * i;
        }
        for(j = 0; j <= rep_period; j++){
            matDP[ 0 * next + j] = INDEL_PENALTY * j;
        }
        // Repeat
        for(i = 0; i < rep_period; i++){
            for(j = 0; j < rep_period; j++ ){
                int gain_or_penalty;
                if( orgInputString[ i + max_pos ] ==
                   orgInputString[ j + sortedString[k] ]){
                    gain_or_penalty = MATCH_GAIN;
                }else{
                    gain_or_penalty = MISMATCH_PENALTY;
                }
                matDP[ (i+1) * next + (j+1)] =
                MAX( matDP[  i    * next + j ]     + gain_or_penalty,
                MAX( matDP[  i    * next + (j+1) ] + INDEL_PENALTY,
                     matDP[ (i+1) * next + j ]     + INDEL_PENALTY ) );
            }
        }
        
#ifdef DEBUG_progressive_multiple_alignment0
        printf("one string = %i\tmax_val = %i\n", sortedString[k],  matDP[rep_period * next + rep_period]);
        for(int i = 0; i <= rep_period; i++){
            for(int j = 0; j <= rep_period; j++ ){
                printf("%i\t", matDP[i * next + j]);
            }
            printf("\n");
        }
#endif
        
        // traceback and revise the consensus alignment
        i = rep_period-1;
        j = rep_period-1;
        int max_val = matDP[rep_period * next + rep_period];
        // max_val is assumed to be matDP[ (i+1) * next + (j+1) ]
        while( i >= 0 ){
            if( i >= 0 && j >= 0 && matDP[ i * next + j ] + MATCH_GAIN == max_val){
                consensus[i][orgInputString[ i + max_pos ]]++;
                max_val -= MATCH_GAIN;
                i--; j--;
            }else if( i >= 0 && j >= 0 && matDP[ i * next + j ] + MISMATCH_PENALTY == max_val){
                consensus[i][orgInputString[ j+sortedString[k] ]]++;
                max_val -= MISMATCH_PENALTY;
                i--; j--;
            }else if( j >= 0 && matDP[ (i+1) * next + j ] + INDEL_PENALTY  == max_val){
                // Memorize it in the gap between i and i+1
                gaps[i+1][orgInputString[ j+sortedString[k] ]]++;
                max_val -= INDEL_PENALTY;
                j--;
            }else if( i >= 0 && matDP[ i * next + (j+1) ] + INDEL_PENALTY == max_val){
                consensus[i][4]++;  // Insert a gap into the reference
                max_val -= INDEL_PENALTY;
                i--;
            }else{
                fprintf(stderr, "fatal error in progressive multiple alignment DP\n");
                exit(EXIT_FAILURE);
            }
        }
        
    }
    // Obtain the consensus string and set it to rep_unit_string with actual_rep_period
    // The current implementation is naive and needs improvement.
    int rep_unit_string[MAX_PERIOD];
    for(int i=0; i < rep_period; i++){
        rep_unit_string[i] = 0;         // default 0=a
        int tmp_max = consensus[i][0];
        
        for(int j=1; j<4; j++){
            if(tmp_max < consensus[i][j]){
                rep_unit_string[i] = j;
                tmp_max = consensus[i][j];
            }
        }
    }
    int actual_rep_period = rep_period;
    
#ifdef DEBUG_progressive_multiple_alignment
    printf("consensus\n");
    for(int i = 0; i < rep_period; i++){
        for(int j = 0; j < 5; j++ ){
            printf("%i\t", consensus[i][j]);
        }
        printf("\n");
    }
    printf("gaps\n");
    for(int i = 0; i < rep_period+1; i++){
        for(int j = 0; j < 5; j++ ){
            printf("%i\t", gaps[i][j]);
        }
        printf("\n");
    }
#endif
    
    rr->actual_rep_period = actual_rep_period;
    rr->ConsensusMethod   = ProgressiveMultipleAlignment;
    if(actual_rep_period > 0)
    {
        wrap_around_DP(rep_unit_string,
                       actual_rep_period,
                       &orgInputString[max_start],
                       (max_end - max_start + 1),
                       rr);
        print_4_decimal_array(rep_unit_string, actual_rep_period, rr->string);
        freq_2mer_array(rep_unit_string, actual_rep_period, rr->freq_2mer);
    }
}

