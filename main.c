//---------------------------------------------------------------------
// Finding tandem repeats in long noizy reads
// Initial codes are developed by Shinichi Morishita and ....
//---------------------------------------------------------------------

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

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
        read_cnt = handle_one_file(inputFile);
        gettimeofday(&e, NULL);
        
#ifdef PRINT_FILE_ETIME
        fprintf(stderr, "The input file name is %s.\n", inputFile);
        fprintf(stderr, "Number of all reads in the input fasta file is %i.\n", read_cnt);
        fprintf(stderr, "time for finding repeats = %lf\n", (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);
#endif
    }
    
    // Remove unqualied reads from repeats_in_all_reads
    int j = 0;
    for(int i=0; i<read_cnt; i++){
        repeat_in_read aRR = repeats_in_all_reads[i];
        float match_ratio;
        if(aRR.actual_repeat_len > 0){
            match_ratio = (float)aRR.Num_matches / aRR.actual_repeat_len;
        }else{
            match_ratio = 0;
        }
        if(MIN_REP_LEN < (aRR.actual_rep_period * aRR.Num_freq_unit)  &&
           MIN_MATCH_RATIO < match_ratio &&
           1 < aRR.Num_freq_unit )
        {
            // Remove unqualied reads from repeats_in_all_reads
            // This overwrite is safe as j <= i.
            repeats_in_all_reads[j] = repeats_in_all_reads[i];
            repeats_in_all_reads[j].ID = j;
            j++;
        }
    }
    read_cnt = j;

    // If it is asked, output a temporary file of reads with repeats
    if(strcmp(argv[1], "-r") == 0){
        for(int i=0; i<read_cnt; i++){
            if(pretty_print == 1){
                pretty_print_one_repeat_in_read(repeats_in_all_reads[i]);
            }else{
                print_one_repeat_in_read(repeats_in_all_reads[i]);
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
        k_means_clustering(read_cnt, pretty_print);
        gettimeofday(&e, NULL);
        
#ifdef PRINT_FILE_ETIME
        fprintf(stderr, "time for clustering repeats = %lf\n", (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);
#endif
    }
    
    //  Free space for global variables in the heap
    free(repeats_in_all_reads);

    return EXIT_SUCCESS;
}
