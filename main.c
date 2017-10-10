//---------------------------------------------------------------------
// Finding tandem repeats in long noizy reads
// Initial codes are developed by Shinichi Morishita and ....
// version 0.6  (September 8, 2017)
//---------------------------------------------------------------------

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

void free_global_variables_and_exit(){
    // If any of the above global variables failed to be allocated in the heap, free other variables and exit.
    if(orgInputString    != NULL){ free(orgInputString); }
    if(inputString       != NULL){ free(inputString); }
    if(count             != NULL){ free(count); }
    if(sortedString      != NULL){ free(sortedString); }
    if(freq_interval_len != NULL){ free(freq_interval_len); }
    if(Kadane_val        != NULL){ free(Kadane_val); }
    if(max_starts        != NULL){ free(max_starts); }
    if(count_period_all  != NULL){ free(count_period_all); }
    if(rep_unit_string   != NULL){ free(rep_unit_string); }
    if(WrapDP            != NULL){ free(WrapDP); }
    fprintf(stderr, "cannot allocate space for one of global variables in the heap.\n");
    exit(EXIT_FAILURE);
}

int handle_one_file(char *inputFile){
    //---------------------------------------------------------------------------
    // Feed a string from a file, convert the string into a series of integers
    //---------------------------------------------------------------------------
    
    FILE *fp = fopen(inputFile, "r");
    if(fp == NULL){
        fprintf(stderr, "fatal error: cannot open %s\n", inputFile);
        exit(EXIT_FAILURE);
    }
    
    int pow4k = 1;
    for(int i=0; i<maxKmer; i++){ pow4k = 4 * pow4k; }  // 4^k
    
    // Allocate the main memory for global variables in the heap
    orgInputString  = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( orgInputString == NULL ){ free_global_variables_and_exit(); }

    inputString     = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( inputString == NULL ){ free_global_variables_and_exit(); }

    count           = (int *)malloc( sizeof(int) * pow4k);
    if( count == NULL ){ free_global_variables_and_exit(); }
    
    sortedString    = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( sortedString == NULL ){ free_global_variables_and_exit(); }
    
    freq_interval_len = (float *)malloc( sizeof(float) * MAX_INPUT_LENGTH);
    if( freq_interval_len == NULL ){ free_global_variables_and_exit(); }
    
    Kadane_val        = (float *)malloc( sizeof(int) * MAX_INPUT_LENGTH);
    if( Kadane_val == NULL ){ free_global_variables_and_exit(); }
    
    max_starts      = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( max_starts == NULL ){ free_global_variables_and_exit(); }
    
    count_period_all= (int *)malloc( sizeof(int) * MAX_PERIOD);
    if( count_period_all == NULL ){ free_global_variables_and_exit(); }
    
    rep_unit_string = (int *)malloc( sizeof(int) * MAX_PERIOD);
    if( rep_unit_string == NULL ){ free_global_variables_and_exit(); }
    
    WrapDP          = (int *)malloc(sizeof(int) * (MAX_PERIOD+1) * (MAX_INPUT_LENGTH+1));
    if( WrapDP == NULL ){ free_global_variables_and_exit(); }
    
    char s[BLK];
    char readID[BLK];
    int i, inputLen, charCode;
    int cnt=0;
    int read_cnt = 0;
    int firstRead = 1;  // 1 means the first read.
    
    while (fgets(s, BLK, fp) != NULL) { // Feed a string of size BLK from fp into string s
        if( MAX_INPUT_LENGTH < cnt){
            fprintf(stderr, "fatal error: the length must be at most %i.\nread ID = %s\n",MAX_INPUT_LENGTH, readID);
            exit(EXIT_FAILURE);
        }
        if(s[0] == '>'){
            if(firstRead == 1){
                firstRead = 0;
            }else{  // Process the previous read if the current read is not the first one.
                inputLen = cnt;
                handle_one_read(readID, inputLen, read_cnt);
                read_cnt++;
            }
            // Feed the header of the read.
            for(i=1; s[i]!='\0' && s[i]!='\n'; i++){ readID[i-1] = s[i]; }
            readID[i-1] = '\0';
            cnt = 0;
        }else{
            for(i=0; s[i]!='\0' && s[i]!='\n'; i++){
                switch(s[i]){
                    case 'A':
                    case 'a':
                        charCode = 0; break;
                    case 'C':
                    case 'c':
                        charCode = 1; break;
                    case 'G':
                    case 'g':
                        charCode = 2; break;
                    case 'T':
                    case 't':
                        charCode = 3; break;
                    default:
                        fprintf(stderr, "Invalid character: %c \n", s[i]); exit(EXIT_FAILURE);
                }
                orgInputString[cnt] = charCode;
                cnt++;  // Count up here.
            }
        }
    }
    inputLen = cnt;
    handle_one_read(readID, inputLen, read_cnt);
    read_cnt++;
    
    fclose(fp);
    
    free(orgInputString);
    free(inputString);
    free(count);
    free(sortedString);
    free(freq_interval_len);
    free(Kadane_val);
    free(max_starts);
    free(count_period_all);
    free(rep_unit_string);
    free(WrapDP);
    
    return(read_cnt);
}

void print_error_message(){
    fprintf(stderr, "Arguments must be of the form -a|-r|-c <file name> (-pp option)\n");
    fprintf(stderr, "-a: Perform All steps. Feed <fasta file>, detect repeats, and print <clustering results>. \n");
    fprintf(stderr, "-r: Find Repeats in reads only. Feed <fasta file> and print <reads with repeats>. \n");
    fprintf(stderr, "-c: Cluster reads with repeats. Feed <reads with repeats> and print <clustering results>. \n");
    fprintf(stderr, "-pp: Pretty print of reads with repeats for better readability and debugging \n");
}

int main(int argc, char *argv[])
{
    // <command> -a <fasta file>    Feed <fasta file>, detect repeats, and print <clustering results>.
    // <command> -r <fasta file>    Feed <fasta file> and print <reads with repeats>.
    // <command> -c <reads with repeats>
    //                              Feed <reads with repeats> and print <clustering results>.
    
    char *inputFile;
    int pretty_print = 0;
    
    if(argc == 3 && (strcmp(argv[1], "-a") == 0 || strcmp(argv[1], "-r") == 0 || strcmp(argv[1], "-c") == 0)){  // Valid arguments
    }else if(argc == 4 && (strcmp(argv[1], "-a") == 0 || strcmp(argv[1], "-r") == 0 || strcmp(argv[1], "-c") == 0) && (strcmp(argv[3], "-pp") == 0)){
            //  Use pretty print if (strcmp(argv[3], "-pp") == 0
        pretty_print = 1;
    }else{
        print_error_message();
        exit(EXIT_FAILURE);
    }
    
    struct timeval s, e;

    //  Allocate space for global variables in the heap
    repeats_in_all_reads = (repeat_in_read*) malloc(MAX_NUM_READS*sizeof(repeat_in_read));
    if(repeats_in_all_reads == NULL){
        fprintf(stderr, "Fatal error: cannot allocate a space for repeats_in_all_reads\n");
        exit(EXIT_FAILURE);
    }
    
    // Process one file to associate reads with tandem repeats
    int read_cnt;
    if(strcmp(argv[1], "-a") == 0 || strcmp(argv[1], "-r") == 0){
        gettimeofday(&s, NULL);
        
        inputFile = argv[2];
        fprintf(stderr, "The input file name is %s\n", inputFile);
        
        read_cnt = handle_one_file(inputFile);
        
        gettimeofday(&e, NULL);
        fprintf(stderr, "time for finding repeats = %lf\n", (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);
    }

    // If it is asked, output a temporary file of reads with repeats
    if(strcmp(argv[1], "-r") == 0){
        for(int i=0; i<read_cnt; i++){
            repeat_in_read aRR = repeats_in_all_reads[i];
            float match_ratio = (float)aRR.Num_matches / aRR.actual_repeat_len;
            
            if(MIN_REP_LEN < (aRR.actual_rep_period * aRR.Num_freq_unit)  &&
               MIN_MATCH_RATIO < match_ratio && 1 < aRR.Num_freq_unit )
            {
                if(pretty_print == 1){
                    pretty_print_one_repeat_in_read(repeats_in_all_reads[i]);
                }else{
                    print_one_repeat_in_read(repeats_in_all_reads[i]);
                }
            }
            if(i % BLK == 0){
                fflush( stdout );
            }
        }
    }
    
    // Cluster reads with repeats into similar groups
    if(strcmp(argv[1], "-c") == 0){
        read_cnt = feed_rr_into_repeats_in_all_reads(argv[2]);
    }
    
    if(strcmp(argv[1], "-a") == 0 || strcmp(argv[1], "-c") == 0){
        gettimeofday(&s, NULL);
        
        int Num_qualified_TRs = k_means_clustering(read_cnt);
        
        for(int j = 0; j < Num_qualified_TRs; j++){
            print_one_TR_with_read(TR_list[j], pretty_print);
            if(j % BLK == 0){
                fflush( stdout );
            }
        }
        
        gettimeofday(&e, NULL);
        fprintf(stderr, "time for clustering repeats = %lf\n", (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);
    }
    
    //  Free space for global variables in the heap
    free(repeats_in_all_reads);

    return EXIT_SUCCESS;
}
