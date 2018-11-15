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
    fprintf(stderr, "Arguments must be of the form -m|-r <file name> \n");
    fprintf(stderr, "-m: Output multiple occurrences of different tandem repeats. \n");
    fprintf(stderr, "-r: Output a single tandem repeat. \n");
}

int main(int argc, char *argv[])
{
    // <command> -m <fasta file>    Feed <fasta file>, detect repeats, and print <clustering results>.
    // <command> -r <fasta file>    Feed <fasta file> and print <reads with repeats>.

    
    char *inputFile;
    int print_multiple_TR;
    
    if(argc == 3 && (strcmp(argv[1], "-m") == 0 || strcmp(argv[1], "-r") == 0)){
        print_multiple_TR = 0;
        if(strcmp(argv[1], "-m") == 0 ){
            print_multiple_TR = 1;
        }
    }else{
        print_error_message();
        exit(EXIT_FAILURE);
    }
    
    // Process one file to associate reads with tandem repeats
    struct timeval s, e;
    int read_cnt;

    gettimeofday(&s, NULL);
    
    inputFile = argv[2];
    fprintf(stderr, "The input file name is %s.\n", inputFile);
             
    read_cnt = handle_one_file(inputFile, print_multiple_TR);
             
    fprintf(stderr, "Number of all reads in the input fasta file %s is %i.\n", inputFile, read_cnt);
    
    gettimeofday(&e, NULL);
    fprintf(stderr, "time for finding repeats = %lf\n", (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);
    
    return EXIT_SUCCESS;
}
