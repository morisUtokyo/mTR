//---------------------------------------------------------------------
// Finding tandem repeats in long noizy reads
// Initial codes are developed by Shinichi Morishita and ....
//---------------------------------------------------------------------

//#include <sys/time.h>
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
    time_all = 0; time_memory = 0; time_range = 0; time_period = 0;
    
    struct timeval s, e;
    gettimeofday(&s, NULL);
    
    inputFile = argv[2];
    fprintf(stderr, "The input file name is %s.\n", inputFile);
    int read_cnt = handle_one_file(inputFile, print_multiple_TR);
    fprintf(stderr, "Number of all reads in the input fasta file\n %s is %i.\n", inputFile, read_cnt);
    
    gettimeofday(&e, NULL);
    time_all = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
#ifdef PRINT_COMP_TIME
    fprintf(stderr, "Computational time\n");
    fprintf(stderr, "%f\tall\n",           time_all);
    
    fprintf(stderr, "%f\tallocating memory\n", time_memory);
    
    fprintf(stderr, "%f\tranges\n",         time_range);
    
    fprintf(stderr, "%f\tComputing periods\n", time_period);
    fprintf(stderr, "\t%f\tpreparation\n",    time_predicted_rep_period_and_max_position);
    fprintf(stderr, "\t%f\tcount table generation\n",   time_count_table);
    fprintf(stderr, "\t%f\tDe Bruijn\n",     time_search_De_Bruijn_graph);
    fprintf(stderr, "\t%f\tprogressive\n",   time_progressive_multiple_alignment);
    fprintf(stderr, "\t%f\twrap around\n",   time_wrap_around_DP);
#endif

    
    return EXIT_SUCCESS;
}
