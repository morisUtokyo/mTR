/*
 Copyright (c) 2019, Shinichi Morishita
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 The views and conclusions contained in the software and documentation are those
 of the authors and should not be interpreted as representing official policies,
 either expressed or implied, of the FreeBSD Project.
 */

//#include <sys/time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

void print_error_message(){
    fprintf(stderr, "mTR [-msp] <file name> \n");
    fprintf(stderr, "-m: Output multiple tandem repeats (default setting). \n");
    fprintf(stderr, "-s: Output the longest tandem repeat. \n");
    fprintf(stderr, "-p: Output the alignment between the input sequence and predicted tandem repeat. \n");
}

int main(int argc, char *argv[])
{
    char *inputFile;
    int print_multiple_TR = 1;
    int print_alignment   = 0;

    if(argc == 2){
        print_multiple_TR = 1;
    }else if(argc == 3){
        char *p;
        p=argv[1];
        if(*p == '-' ){
            for(p++; *p != '\0'; p++){
                switch(*p){
                    case 's':   print_multiple_TR = 0;
                                MIN_MAX_DI = MIN_MAX_DI_SINGLE;   // setting 0 does not improve the accuracy
                                break;
                    case 'm':   print_multiple_TR = 1;  break;
                                MIN_MAX_DI = 0;
                    case 'p':   print_alignment   = 1;  break;
                    default:    MIN_MAX_DI = 0; // multiple mode
                }
            }
        }
    }else{
        print_error_message();
        exit(EXIT_FAILURE);
    }
    
    // Process one file to associate reads with tandem repeats
    time_all = 0; time_memory = 0; time_range = 0; time_period = 0;
    time_initialize_input_string = 0;
    time_count_table = 0; time_search_De_Bruijn_graph = 0;
    time_wrap_around_DP = 0; time_chaining = 0;
    query_counter = 0;
    
    
    struct timeval s, e;
    gettimeofday(&s, NULL);
    
    inputFile = argv[2];
    fprintf(stderr, "The input file name is %s.\n", inputFile);
    int read_cnt = handle_one_file(inputFile, print_multiple_TR, print_alignment);
    fprintf(stderr, "Number of reads is %i.\n", read_cnt);
    
    gettimeofday(&e, NULL);
    time_all = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
#ifdef PRINT_COMP_TIME
    fprintf(stderr, "Computational time\n");
    fprintf(stderr, "%f\tall\n",           time_all);
    fprintf(stderr, "%f\tallocating memory\n", time_memory);
    fprintf(stderr, "%f\tranges\n",         time_range);
    
    fprintf(stderr, "%f\tComputing periods\n", time_period);
    fprintf(stderr, "\t%f\tInitialize the input\n", time_initialize_input_string);
    fprintf(stderr, "\t%f\tcount table generation\n",   time_count_table);
    fprintf(stderr, "\t%f\tDe Bruijn\n",     time_search_De_Bruijn_graph);
    fprintf(stderr, "\t%f\twrap around\n",   time_wrap_around_DP);
    fprintf(stderr, "\t%f\tchaining\n",   time_chaining);
    fprintf(stderr, "\t%i\tCount of queries\n", query_counter);
#endif

    
    return EXIT_SUCCESS;
}
