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
    fprintf(stderr, "Arguments must be of the form (-m|-s) <file name> \n");
    fprintf(stderr, "-m: Output multiple tandem repeats (default setting). \n");
    fprintf(stderr, "-s: Output the longest tandem repeat. \n");
}

int main(int argc, char *argv[])
{
    // <command> -m <fasta file>    Feed <fasta file>, detect repeats, and print <clustering results>.
    // <command> -s <fasta file>    Feed <fasta file> and print <reads with repeats>.

    
    char *inputFile;
    int print_multiple_TR = 1;

    if(argc == 2){
        print_multiple_TR = 1;
    }else if(argc == 3){
        char *p;
        p=argv[1];
        if(*p == '-' ){
            for(p++; *p != '\0'; p++){
                switch(*p){
                    case 's':   print_multiple_TR = 0; break;
                }
            }
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
    fprintf(stderr, "Number of reads is %i.\n", read_cnt);
    
    gettimeofday(&e, NULL);
    time_all = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
#ifdef PRINT_COMP_TIME
    fprintf(stderr, "Computational time\n");
    fprintf(stderr, "%f\tall\n",           time_all);
    
    fprintf(stderr, "%f\tallocating memory\n", time_memory);
    
    fprintf(stderr, "%f\tranges\n",         time_range);
    
    fprintf(stderr, "%f\tComputing periods\n", time_period);
    fprintf(stderr, "\t%f\tcount table generation\n",   time_count_table);
    fprintf(stderr, "\t%f\tDe Bruijn\n",     time_search_De_Bruijn_graph);
    fprintf(stderr, "\t%f\twrap around\n",   time_wrap_around_DP);
#endif

    
    return EXIT_SUCCESS;
}
