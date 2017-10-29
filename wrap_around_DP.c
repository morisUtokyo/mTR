//
//  wrap_around_DP.c
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

// We can select either local or global alignment by setting mode to Local_alignment or Global_alignment.
// We are using local alignment as default.
void wrap_around_DP_sub(int *rep_unit, int unit_len, int *rep, int rep_len, repeat_in_read *rr, int mode){
    
    int i, j;
    int next = unit_len+1;
    
#ifdef DEBUG_algorithm_wrap_around
    printf("unit_len = %i, \t rep_len = %i\n", unit_len, rep_len);
#endif
    
    // Initialization
    for(j=0; j<=rep_len; j++){   // Scan rep_unit
        if(mode == Global_alignment){
            // global alignment
            WrapDP[next*0 + j] = INDEL_PENALTY * j;
        }else{
            // local alignment
            WrapDP[next*0 + j] = 0;
        }
    }
    int max_wrd = 0;
    int max_i = 0;
    int max_j = 0;
    for(i=1; i<=rep_len; i++){        // Scan rep
        for(j=1; j<=unit_len; j++){   // Scan rep_unit
            if(rep[i-1] == rep_unit[j-1]){    // 0-origin index !!!!
                WrapDP[next*i + j] =
                WrapDP[next*(i-1) + j-1] + MATCH_GAIN;
            }else{
                if(mode == Global_alignment){
                    // global alignment
                    if(j == 1){
                        WrapDP[next*i + j] = MAX(
                                  WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                  WrapDP[next*(i-1) + j]   + INDEL_PENALTY);
                    }else{
                        WrapDP[next*i + j] = MAX( MAX(
                                  WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY,
                                  WrapDP[next*i + j-1]     + INDEL_PENALTY),
                                  WrapDP[next*(i-1) + j]   + INDEL_PENALTY);
                    }
                }else{
                    // local alignment
                    if(j == 1){
                        WrapDP[next*i + j] = MAX( MAX( 0,
                                  WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY),
                                  WrapDP[next*(i-1) + j]   + INDEL_PENALTY);
                    }else{
                        WrapDP[next*i + j] = MAX( MAX( MAX( 0,
                                  WrapDP[next*(i-1) + j-1] + MISMATCH_PENALTY),
                                  WrapDP[next*i + j-1]     + INDEL_PENALTY),
                                  WrapDP[next*(i-1) + j]   + INDEL_PENALTY);
                    }
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
    
    // trace back
    int Num_matches = 0;
    int Num_mismatches = 0;
    int Num_insertions = 0;
    int Num_deletions  = 0;
    int Num_scanned_unit = 0;
    
    i = max_i;
    j = max_j;
    if(j == 0){ j = unit_len; }
    
    while(i > 0){
        if( mode == Local_alignment && max_wrd <= 0){ break; }
        // If we use local alignment, exit here when max_wrd is equal to or less than 0.
        
        if( i > 0 && j > 0 && max_wrd == WrapDP[next*(i-1) + j-1]+ MATCH_GAIN){
            max_wrd -= MATCH_GAIN;
            i--; j--;
            Num_matches++;
            Num_scanned_unit++;
        }else if( i > 0 && j > 0 && max_wrd == WrapDP[next*(i-1) + j-1]+ MISMATCH_PENALTY){
            max_wrd -= MISMATCH_PENALTY;
            i--; j--;
            Num_mismatches++;
            Num_scanned_unit++;
        }else if( i >= 0 && j > 0 && max_wrd == WrapDP[next*i + j-1]    + INDEL_PENALTY){
            max_wrd -= INDEL_PENALTY;
            j--;
            Num_insertions++;
            Num_scanned_unit++;
        }else if( i > 0 && j >= 0 && max_wrd == WrapDP[next*(i-1) + j] + INDEL_PENALTY){
            max_wrd -= INDEL_PENALTY;
            i--;
            Num_deletions++;
        }else{
            fprintf(stderr, "fatal error in wrap-around DP\n");
            exit(EXIT_FAILURE);
        }
        if(j == 0){
            j = unit_len;
        }
    }
    int actual_repeat_len   = max_i - i;
    rr->actual_repeat_len   = actual_repeat_len;
    rr->Num_freq_unit       = (int)Num_scanned_unit/unit_len;
    rr->Num_matches         = Num_matches;
    rr->Num_mismatches      = Num_mismatches;
    rr->Num_insertions      = Num_insertions;
    rr->Num_deletions       = Num_deletions;
    
#ifdef DEBUG_algorithm_wrap_around
    printf("Repeat length =%i\t Matches = %1.3f\t Mismatches = %1.3f\t Insertions = %1.3f\t Deletions = %1.3f\t max value = %i\n",
           actual_repeat_len,
           (float)Num_matches/actual_repeat_len,
           (float)Num_mismatches/actual_repeat_len,
           (float)Num_insertions/actual_repeat_len,
           (float)Num_deletions/actual_repeat_len,
           Num_matches*MATCH_GAIN + Num_mismatches*MISMATCH_PENALTY + Num_insertions * INDEL_PENALTY + Num_deletions * INDEL_PENALTY);
#endif
    
}

void wrap_around_DP(int *rep_unit, int unit_len, int *rep, int rep_len, repeat_in_read *rr){

    int mode = Wrap_around_DP_mode;   // Set to either global_alignment or local_alignment in mTR.h
    wrap_around_DP_sub(rep_unit, unit_len, rep, rep_len, rr, mode);
    
}
