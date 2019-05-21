//
//  wrap_around_DP.c
//  
//
//  Created by Shinichi Morishita on 2017/10/06.
//

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

char base_int2char(int val){
    switch(val){
        case 0: return('A');
        case 1: return('C');
        case 2: return('G');
        case 3: return('T');
        default: fprintf(stderr, "fatal input char %i\n", val); exit(EXIT_FAILURE);
    }
}

int base_char2int(char c){
    switch(c){
        case 'A': return(0);
        case 'C': return(1);
        case 'G': return(2);
        case 'T': return(3);
        default: fprintf(stderr, "fatal input char %c", c); exit(EXIT_FAILURE);
    }
}

void pretty_print_alignment(int *rep_unit, int unit_len, int rep_start, int rep_end){
    
    int *rep;
    rep = &orgInputString[rep_start];
    int rep_len = rep_end - rep_start + 1;
    
    int i, j;
    int next = unit_len+1;
    
    // Initialization
    for(j=0; j<=unit_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = 0; // local alignment
    }
    int max_wrd = 0;
    int max_i = 0;
    int max_j = 0;
    int val_match, val_mismatch, val_insertion, val_deletion;
    for(i=1; i <= rep_len; i++){        // Scan rep
        for(j=1; j<=unit_len; j++){   // Scan rep_unit
            if( WrapDPsize <= next*i + j ){
                fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
                exit(EXIT_FAILURE);
            }
            val_match       = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
            val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
            val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
            
            if(rep[i] == rep_unit[j]){    // *1*-origin index !!!!
                WrapDP[next*i + j] = val_match;
            }else{
                if(j > 1){
                    val_deletion = WrapDP[next*i + j-1] - INDEL_PENALTY;
                    WrapDP[next*i + j] =
                    MAX(0, MAX( MAX( val_mismatch, val_insertion), val_deletion));
                }else{
                    WrapDP[next*i + j] =
                    MAX(0, MAX( val_mismatch, val_insertion));
                }
            }
            if(max_wrd < WrapDP[next*i+j]){
                max_wrd = WrapDP[next*i+j];
                max_i = i;
                max_j = j;
            }
        }
        // wrap around
        WrapDP[next*i + 0] = WrapDP[next*i + unit_len];
    }
    
    // trace back the optimal alignment while storing it in the data structure "alignment"
    
    int Num_matches = 0;
    int Num_mismatches = 0;
    int Num_insertions = 0;
    int Num_deletions  = 0;
    int pos = 0;
    
    i = max_i;
    j = max_j;
    if(j == 0){ j = unit_len; } // 1-origin index
    
    while(i > 0 && WrapDP[next*i + j] > 0){                 // global alignment
        val_match       = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
        val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
        val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
        val_deletion    = WrapDP[next*i + j-1]      - INDEL_PENALTY;
        
        if( max_wrd == val_match && rep[i] == rep_unit[j]){
            alignment_input[pos]   = base_int2char(rep[i]);
            alignment_symbols[pos] = '|';
            alignment_repeats[pos] = base_int2char(rep_unit[j]);
            
            max_wrd -= MATCH_GAIN;
            i--; j--;
            Num_matches++;
            pos++;
        }else if( max_wrd == val_mismatch  && rep[i] != rep_unit[j]){     // mismatch
            alignment_input[pos]   = base_int2char(rep[i]);
            alignment_symbols[pos] = ' ';
            alignment_repeats[pos] = base_int2char(rep_unit[j]);
            
            max_wrd += MISMATCH_PENALTY;
            i--; j--;
            Num_mismatches++;
            pos++;
        }else if( max_wrd == val_deletion){     // deletion
            alignment_input[pos]   = '-';
            alignment_symbols[pos] = ' ';
            alignment_repeats[pos] = base_int2char(rep_unit[j]);

            max_wrd += INDEL_PENALTY;
            j--;
            Num_deletions++;    // Num_insertions++;
            pos++;
        }else if( max_wrd == val_insertion){    // insertion
            alignment_input[pos]   = base_int2char(rep[i]);
            alignment_symbols[pos] = ' ';
            alignment_repeats[pos] = '-';

            max_wrd += INDEL_PENALTY;
            i--;
            Num_insertions++;   // Num_deletions++;
            pos++;
        }else if( max_wrd == 0){
            break;
        }else{
             fprintf(stderr, "fatal error in wrap-around DP max_wrd = %i\n", max_wrd);
             exit(EXIT_FAILURE);
        }
        if(j == 0){
            j = unit_len;
        }
    }
    
    //printf("rep_start = %i\trep_end = %i\t#Matches=%i\t#Mismatches=%i\t#Insertions=%i\t#Deletions=%i\n\n", rep_start, rep_end, Num_matches, Num_mismatches, Num_insertions, Num_deletions);
    

    // print backward from the end of alignment arrays
    for(int i_start = pos-1; 0 <= i_start; i_start -= ALIGNMENT_WIDTH_PRINTING){
        int i_end;
        if(-1 <= i_start - ALIGNMENT_WIDTH_PRINTING){
            i_end = i_start - ALIGNMENT_WIDTH_PRINTING;
        }else{
            i_end = -1;
        }
        for(int i = i_start; i_end < i; i--){
            printf("%c", alignment_input[i]);
        }
        printf("\n");
        //printf("\t%i-%i\n", rep_start + (pos - 1 - i_start), rep_start + (pos - 1 - i_end - 1));
        
        for(int i = i_start; i_end < i; i--){
            printf("%c", alignment_symbols[i]);
        }
        printf("\n");
        for(int i = i_start; i_end < i; i--){
            printf("%c", alignment_repeats[i]);
        }
        printf("\n\n");
    }
}



// The repeat unit is the input sequence.
// The number of deletions is represented as a positive integer.

// We used a local alignment algorithm in place of a global alignment that outputs lengthy alignments.
void wrap_around_DP(int *rep_unit,          int unit_len,
                    int query_start,        int query_end,
                    int *actual_start,      int *actual_end,
                    int *return_rep_len,    int *return_freq_unit,
                    int *return_matches,    int *return_mismatches,
                    int *return_insertions, int *return_deletions){
    
    struct timeval s, e;
    gettimeofday(&s, NULL);

    // Initialization
    int *rep;
    rep = &orgInputString[query_start];
    int rep_len = query_end - query_start + 1;
    
    int i, j;
    int next = unit_len+1;
    
    for(j=0; j<=unit_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = 0; // local alignment
    }
    int max_wrd = 0;
    int max_i = 0;
    int max_j = 0;
    int val_match, val_mismatch, val_insertion, val_deletion;
    for(i=1; i <= rep_len; i++){        // Scan rep
        for(j=1; j<=unit_len; j++){   // Scan rep_unit
            if( WrapDPsize <= next*i + j ){
                fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
                exit(EXIT_FAILURE);
            }
            val_match       = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
            val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
            val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
            
            if(rep[i] == rep_unit[j]){    // *1*-origin index !!!!
                WrapDP[next*i + j] = val_match;
            }else{
                if(j > 1){
                    val_deletion = WrapDP[next*i + j-1] - INDEL_PENALTY;
                    WrapDP[next*i + j] =
                        MAX(0, MAX( MAX( val_mismatch, val_insertion), val_deletion));
                }else{
                    WrapDP[next*i + j] =
                        MAX(0, MAX( val_mismatch, val_insertion));
                }
            }
            if(max_wrd < WrapDP[next*i+j]){
                max_wrd = WrapDP[next*i+j];
                max_i = i;
                max_j = j;
            }
        }
        // wrap around
        WrapDP[next*i + 0] = WrapDP[next*i + unit_len];
    }
    
    // trace back the optimal alignment while storing it in the data structure "alignment"
    int Num_matches = 0;
    int Num_mismatches = 0;
    int Num_insertions = 0;
    int Num_deletions  = 0;
    int Num_scanned_unit = 0;
    
    i = max_i;
    j = max_j;
    if(j == 0){ j = unit_len; } // 1-origin index
    
    while(i > 0 && WrapDP[next*i + j] > 0){                 // global alignment
        val_match       = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
        val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
        val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
        val_deletion    = WrapDP[next*i + j-1]      - INDEL_PENALTY;
        
        if( max_wrd == val_match          && rep[i] == rep_unit[j]){
            max_wrd -= MATCH_GAIN;
            i--; j--;
            Num_matches++;
            Num_scanned_unit++;
        }else if( max_wrd == val_mismatch && rep[i] != rep_unit[j]){     // mismatch
            max_wrd += MISMATCH_PENALTY;
            i--; j--;
            Num_mismatches++;
            Num_scanned_unit++;
        }else if( max_wrd == val_deletion){     // deletion
            max_wrd += INDEL_PENALTY;
            j--;
            Num_deletions++;    // Num_insertions++;
            Num_scanned_unit++;
        }else if( max_wrd == val_insertion){    // insertion
            max_wrd += INDEL_PENALTY;
            i--;
            Num_insertions++;
            //Num_scanned_unit++;       // The base of the repeat unit is skipped.
        }else if( max_wrd == 0){
            break;
        }else{
            fprintf(stderr, "fatal error in wrap-around DP max_wrd = %i\n", max_wrd);
            exit(EXIT_FAILURE);
        }
        if(j == 0){
            j = unit_len;
        }
    }
 
    int actual_repeat_len = max_i - i + 1;
    
    *actual_start       = query_start + i;
    *actual_end         = query_start + max_i;
    *return_rep_len     = actual_repeat_len;
    *return_freq_unit   = (int)Num_scanned_unit/unit_len;
    *return_matches     = Num_matches;
    *return_mismatches  = Num_mismatches;
    *return_insertions  = Num_insertions;
    *return_deletions   = Num_deletions;
    
 
    gettimeofday(&e, NULL);
    time_wrap_around_DP += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
}

