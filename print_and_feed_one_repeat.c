//
//  print_one_repeat.c
//  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

int feed_rr_into_repeats_in_all_reads(char *inputFile){
    FILE *fp = fopen(inputFile, "r");
    if(fp == NULL){
        fprintf(stderr, "fatal error: cannot open %s\n", inputFile);
        exit(EXIT_FAILURE);
    }
    char s[BLK];
    int i = 0;
    repeat_in_read tmp;
    
    while( fscanf(fp, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
             &tmp.ID,
             tmp.readID,
             &tmp.inputLen,
             &tmp.max_start,
             &tmp.max_end,
             &tmp.actual_repeat_len,
             &tmp.rep_period,
             &tmp.actual_rep_period,
             &tmp.Num_freq_unit,
             &tmp.Num_matches,
             &tmp.Num_mismatches,
             &tmp.Num_insertions,
             &tmp.Num_deletions,
             &tmp.Kmer,
             &tmp.ConsensusMethod,
             tmp.string,
             &tmp.freq_2mer[0], &tmp.freq_2mer[1],
             &tmp.freq_2mer[2], &tmp.freq_2mer[3],
             &tmp.freq_2mer[4], &tmp.freq_2mer[5],
             &tmp.freq_2mer[6], &tmp.freq_2mer[7],
             &tmp.freq_2mer[8], &tmp.freq_2mer[9],
             &tmp.freq_2mer[10],&tmp.freq_2mer[11],
             &tmp.freq_2mer[12],&tmp.freq_2mer[13],
             &tmp.freq_2mer[14],&tmp.freq_2mer[15]
                  ) != EOF )
    {
        repeats_in_all_reads[i] = tmp;
        repeats_in_all_reads[i].ID = i;
        i++;
    }
    
    fclose(fp);
    return(i);
}

void print_one_repeat_in_read(repeat_in_read rr){
    char message[2000] = "";
    sprintf(message, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s",
            rr.ID,
            rr.readID,
            rr.inputLen,
            rr.max_start,
            rr.max_end,
            rr.actual_repeat_len,
            rr.rep_period,
            rr.actual_rep_period,
            rr.Num_freq_unit,
            rr.Num_matches,
            rr.Num_mismatches,
            rr.Num_insertions,
            rr.Num_deletions,
            rr.Kmer,
            rr.ConsensusMethod,
            rr.string
            );
    char tmp_msg[6]="";
    for(int i=0; i<16; i++){
        sprintf(tmp_msg, "\t%d", rr.freq_2mer[i]);
        strcat(message, tmp_msg);
    }
    printf("%s\n", message);
}


void pretty_print_one_repeat_in_read(repeat_in_read rr){
    
    char message[2000] = "";
    char consensus[100];
    if(rr.ConsensusMethod == 0){
        strcpy(consensus, "Progressive MA ");
    }else{
        strcpy(consensus, "De Bruijn graph");
    }
    
    if(rr.actual_rep_period > 0){
        sprintf(message, "ID = %i,\treadID = %s,\tREAD_length = %i,\tREPEAT_start = %i,\tend = %i,\tlen(act) = %i,\tUNIT_len(est) = %i,\tlen(act) = %i,\t#Units = %i,\tMatches = %i,\tMismacthes = %i,\tIns = %i,\tDel = %i,\tKmer = %i,\tMethod = %s,\tUnit_string= %s,\t2mer_freq = ",
                rr.ID,
                rr.readID,
                rr.inputLen,
                rr.max_start,
                rr.max_end,
                rr.actual_repeat_len,
                rr.rep_period,
                rr.actual_rep_period,
                rr.Num_freq_unit,
                rr.Num_matches,
                rr.Num_mismatches,
                rr.Num_insertions,
                rr.Num_deletions,
                rr.Kmer,
                consensus,
                rr.string
                );
        char tmp_msg[6]="";
        for(int i=0; i<16; i++){
            sprintf(tmp_msg, "%i,", rr.freq_2mer[i]);
            strcat(message, tmp_msg);
        }
    }else{
        sprintf(message, "ID = %i,\treadID = %s,\tREAD_length = %i,\tREPEAT_start = %i,\tend = %i,\tlen(est) = %i,\tlen(act) = NA,\t\tUNIT len(est) = %i,\tlen(act) = NA for Kmer = %i",
                rr.ID,
                rr.readID,
                rr.inputLen,
                rr.max_start,
                rr.max_end,
                (rr.max_end - rr.max_start + 1),
                rr.rep_period,
                rr.Kmer
                );
    }
    printf("%s\n", message);
}

