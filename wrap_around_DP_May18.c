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

void print_one_base(int val){
    switch(val){
        case 0: printf("A"); break;
        case 1: printf("C"); break;
        case 2: printf("G"); break;
        case 3: printf("T"); break;
        default: fprintf(stderr, "fatal input char %i", val%4); exit(EXIT_FAILURE);
    }
}

void print_alignment(int *rep_unit, int unit_len, int *rep, int rep_len){
    
    struct timeval s, e;
    gettimeofday(&s, NULL);
    
    int i, j;
    int next = unit_len+1;
    
    
    // Initialization
#ifdef LOCAL_ALIGNMENT
    for(j=0; j<=rep_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = 0; // local alignment
    }
#else
    for(j=0; j<=rep_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = INDEL_PENALTY * j; // global alignment
    }
#endif
    
    int max_wrd = 0;
    int max_i = 0;
    int max_j = 0;
    for(i=1; i<=rep_len; i++){        // Scan rep
        for(j=1; j<=unit_len; j++){   // Scan rep_unit
            if( WrapDPsize <= next*i + j ){
                fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
                exit(EXIT_FAILURE);
            }
            
            if(rep[i-1] == rep_unit[j-1]){    // 0-origin index !!!!
                WrapDP[next*i + j] =
                WrapDP[next*(i-1) + j-1] + MATCH_GAIN;
            }else{
                if(j == 1){
#ifdef LOCAL_ALIGNMENT
                    WrapDP[next*i + j] = MAX(0, MAX(
                                                    WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                                    WrapDP[next*(i-1) + j]   + INDEL_PENALTY));
#else
                    WrapDP[next*i + j] = MAX(
                                             WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                             WrapDP[next*(i-1) + j]   + INDEL_PENALTY);
#endif
                }else{
#ifdef LOCAL_ALIGNMENT
                    WrapDP[next*i + j] = MAX(0, MAX( MAX(
                                                    WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                                    WrapDP[next*i + j-1]     + INDEL_PENALTY),
                                                    WrapDP[next*(i-1) + j]   + INDEL_PENALTY));
#else
                    WrapDP[next*i + j] = MAX( MAX(
                                                    WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                                    WrapDP[next*i + j-1]     + INDEL_PENALTY),
                                                    WrapDP[next*(i-1) + j]   + INDEL_PENALTY);
#endif
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
        if( i > 0 && j > 0 && max_wrd == WrapDP[next*(i-1) + j-1] + MATCH_GAIN){
            max_wrd -= MATCH_GAIN;
            alignment[i-1] = j-1;
            i--; j--;
            Num_matches++;
            Num_scanned_unit++;
        }else if( i > 0 && j > 0 && max_wrd == WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY){
            max_wrd -= MISMATCH_PENALTY;
            alignment[i-1] = j-1;
            i--; j--;
            Num_mismatches++;
            Num_scanned_unit++;
        }else if( i >= 0 && j > 0 && max_wrd == WrapDP[next*i + j-1] + INDEL_PENALTY){
            // deletion
            max_wrd -= INDEL_PENALTY;
            j--;
            Num_deletions++;    // Num_insertions++;
            Num_scanned_unit++;
        }else if( i > 0 && j >= 0 && max_wrd == WrapDP[next*(i-1) + j] + INDEL_PENALTY){
            // insertion
            max_wrd -= INDEL_PENALTY;
            alignment[i-1] = j;
            i--;
            Num_insertions++;   // Num_deletions++;
        }else{
            fprintf(stderr, "fatal error in wrap-around DP\n");
            exit(EXIT_FAILURE);
        }
        if(j == 0){
            j = unit_len;
        }
    }
    int actual_start       = i+1;
    int actual_end         = max_i;
    
    int diff, i_end;
    int print_width = 50;
    
    printf("%i\n", unit_len);
    for(int i_start = actual_start; i_start < actual_end; i_start += print_width){
        if(i_start + print_width < actual_end){
            i_end = i_start + print_width;
        }else{
            i_end = actual_end;
        }
        // print input sequence
        for(int i = i_start; i < i_end; i++){
             diff = (alignment[i]  - alignment[i-1] + unit_len) % unit_len;
            if(diff == 1 || diff == 0){
                // one match, mismatch or insertion
                print_one_base(rep[i-1]);
            }else if(diff > 1){
                // deletions
                for(int j=0; j<diff-1; j++){
                    printf("-");
                }
            }else{
                fprintf(stderr, "fatal error in the alignment data structure: diff = %i", diff);
                exit(EXIT_FAILURE);
            }
        }
        printf("\n");
        // print symbols of matches, mismatches, insertions, and deletions
        for(int i = i_start; i < i_end; i++){
             diff = (alignment[i]  - alignment[i-1] + unit_len) % unit_len;
            if(diff == 1){
                // one match or mismatch
                if(rep[i-1] == rep_unit[alignment[i-1]]){
                    printf("|");
                }else{
                    printf(" ");
                }
            }else if(diff > 1){
                // deletions
                for(int j=0; j<diff-1; j++){
                    printf(" ");
                }
            }else if(diff == 0){
                // insertions
                printf(" ");
            }else{
                fprintf(stderr, "fatal error in the alignment data structure: diff = %i", diff);
                exit(EXIT_FAILURE);
            }
        }
        printf("\n");
        // print repeat units
        for(int i = i_start; i < i_end; i++){
            diff = (alignment[i]  - alignment[i-1] + unit_len) % unit_len;
            if(diff == 1){
                // one match or mismatch
                print_one_base(rep_unit[alignment[i-1]]);
            }else if(diff > 1){
                // deletions
                for(int j=0; j<diff-1; j++){
                    print_one_base(rep_unit[alignment[i-1+j]]);
                }
            }else if(diff == 0){
                // insertions
                printf("-");
            }else{
                fprintf(stderr, "fatal error in the alignment data structure: diff = %i", diff);
                exit(EXIT_FAILURE);
            }
        }
        printf("\n\n");
    }
}



// The repeat unit is the input sequence.
// The number of deletions is represented as a positive integer.

// We used a local alignment algorithm in place of a global alignment that outputs lengthy alignments.
void wrap_around_DP(int *rep_unit, int unit_len, int *rep, int rep_len,
                    int *actual_start,      int *actual_end,
                    int *return_rep_len,    int *return_freq_unit,
                    int *return_matches,    int *return_mismatches,
                    int *return_insertions, int *return_deletions){
    
    struct timeval s, e;
    gettimeofday(&s, NULL);
    
    int i, j;
    int next = unit_len+1;
    
    
    
#ifdef DEBUG_algorithm_wrap_around
    printf("unit_len = %i, \t rep_len = %i\n", unit_len, rep_len);
#endif
    
    // Initialization
#ifdef LOCAL_ALIGNMENT
    for(j=0; j<=rep_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = 0; // local alignment
    }
#else
    for(j=0; j<=rep_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = INDEL_PENALTY * j; // global alignment
    }
#endif
    
    int max_wrd = 0;
    int max_i = 0;
    int max_j = 0;
    for(i=1; i<=rep_len; i++){        // Scan rep
        for(j=1; j<=unit_len; j++){   // Scan rep_unit
            if( WrapDPsize <= next*i + j ){
                fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
                exit(EXIT_FAILURE);
            }
            
            if(rep[i-1] == rep_unit[j-1]){    // 0-origin index !!!!
                WrapDP[next*i + j] =
                WrapDP[next*(i-1) + j-1] + MATCH_GAIN;
            }else{
                if(j == 1){
#ifdef LOCAL_ALIGNMENT
                    WrapDP[next*i + j] = MAX(0, MAX(
                                             WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                             WrapDP[next*(i-1) + j]   + INDEL_PENALTY));
#else
                    WrapDP[next*i + j] = MAX(
                                         WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                         WrapDP[next*(i-1) + j]   + INDEL_PENALTY);
#endif
                }else{
#ifdef LOCAL_ALIGNMENT
                    WrapDP[next*i + j] = MAX(0, MAX( MAX(
                                          WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                          WrapDP[next*i + j-1]     + INDEL_PENALTY),
                                          WrapDP[next*(i-1) + j]   + INDEL_PENALTY));
#else
                    WrapDP[next*i + j] = MAX( MAX(
                                          WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                          WrapDP[next*i + j-1]     + INDEL_PENALTY),
                                          WrapDP[next*(i-1) + j]   + INDEL_PENALTY);
#endif
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
    
#ifdef DEBUG_algorithm_wrap_around_all
    printf("repeat length = %i\t max value = %i\t max i = %i\t max j =%i\n", rep_len, max_wrd, max_i, max_j);
    
    for(i=0; i<=rep_len; i++){        // Scan rep
        printf("%i\t", i);
        for(j=0; j<=unit_len; j++){   // Scan rep_unit
            printf("%i\t", WrapDP[next*i + j]);
        }
        printf("\n");
    }
#endif
    
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
        if( i > 0 && j > 0 && max_wrd == WrapDP[next*(i-1) + j-1] + MATCH_GAIN){
            max_wrd -= MATCH_GAIN;
            alignment[i] = j;
            i--; j--;
            Num_matches++;
            Num_scanned_unit++;
        }else if( i > 0 && j > 0 && max_wrd == WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY){
            max_wrd -= MISMATCH_PENALTY;
            alignment[i] = j;
            i--; j--;
            Num_mismatches++;
            Num_scanned_unit++;
        }else if( i >= 0 && j > 0 && max_wrd == WrapDP[next*i + j-1] + INDEL_PENALTY){
            // deletion
            max_wrd -= INDEL_PENALTY;
            j--;
            Num_deletions++;    // Num_insertions++;
            Num_scanned_unit++;
        }else if( i > 0 && j >= 0 && max_wrd == WrapDP[next*(i-1) + j] + INDEL_PENALTY){
            // insertion
            max_wrd -= INDEL_PENALTY;
            alignment[i] = j;
            i--;
            Num_insertions++;   // Num_deletions++;
        }else{
            fprintf(stderr, "fatal error in wrap-around DP\n");
            exit(EXIT_FAILURE);
        }
        if(j == 0){
            j = unit_len;
        }
    }
    *actual_start       = i;
    *actual_end         = max_i;
    
    int actual_repeat_len = max_i - i + 1;
    *return_rep_len     = actual_repeat_len;
    *return_freq_unit   = (int)Num_scanned_unit/unit_len;
    *return_matches     = Num_matches;
    *return_mismatches  = Num_mismatches;
    *return_insertions  = Num_insertions;
    *return_deletions   = Num_deletions;
    
#ifdef PRINT_ALIGNMENT
    printf("Repeat length =%i\t Matches = %i\t Mismatches = %i\t Insertions = %i\t Deletions = %i\n",
           actual_repeat_len,
           Num_matches,
           Num_mismatches,
           Num_insertions,
           Num_deletions);
    
    print_alignment(rep_unit, unit_len, rep, rep_len);
#endif
 
    gettimeofday(&e, NULL);
    time_wrap_around_DP += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
}

