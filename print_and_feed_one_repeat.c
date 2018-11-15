//
//  print_one_repeat.c
//
//
//  Created by Shinichi Morishita
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

void print_one_repeat_in_read(repeat_in_read rr){
    printf(
            "%.10s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%s\n",
            //rr.ID,
            rr.readID,
            rr.inputLen,
            rr.max_start,
            rr.max_end,
            rr.actual_repeat_len,
            rr.actual_rep_period,
            rr.Num_freq_unit,
            rr.Num_matches,
            (float)rr.Num_matches/rr.actual_repeat_len,
            rr.Num_mismatches,
            rr.Num_insertions,
            rr.Num_deletions,
            rr.Kmer,
            rr.ConsensusMethod,
            rr.string
            );
    fflush(stdout);
}

