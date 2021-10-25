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
#include <unistd.h>
#include "mTR.h"
#include <getopt.h>

#include <mpi.h>

void print_error_message(){
    fprintf(stderr, "mTR [-acp] [-m ratio] <fasta file name> \n");
    fprintf(stderr, "-a: Output the alignment between the input sequence and predicted tandem repeat. \n");
    fprintf(stderr, "-c: Print the computation time of each step.\n");
    fprintf(stderr, "-m ratio: Give a minimum match ratio ranging from 0 to 1.\n");
    fprintf(stderr, "-p: Use Pearson's correlation coefficient distance in place of Manhattan distance.\n");
}

int main(int argc, char *argv[])
{
    char *inputFile;

    // MPI variables
    int myid, num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // default parameters
    int print_computation_time = 0;
    min_match_ratio = MIN_MATCH_RATIO;
    int print_alignment = 0;
    Manhattan_Distance = 1;
    
    int opt;
    while ((opt = getopt(argc, argv, "acm:p")) != -1) {
        switch(opt){
            case 'a':
                print_alignment = 1;
                break;
            case 'c':
                print_computation_time = 1;
                break;
            case 'm':
                min_match_ratio = atof(optarg);
                //fprintf(stderr, "The input minimum match ratio is %f.\n", min_match_ratio);
                if(0 <= min_match_ratio && min_match_ratio <=1){
                    break;
                }else{
                    fprintf(stderr, "The input minimum match ratio must range from 0 to 1.\n", min_match_ratio);
                    exit(EXIT_FAILURE);
                }
            case 'p':
                Manhattan_Distance = 0;
                fprintf(stderr, "Pearson's correlation coefficient distance in place of Manhattan distance.\n");
                break;
            default:
                print_error_message();
                exit(EXIT_FAILURE);
        }
    }
    if (optind >= argc) {
        fprintf(stderr, "The input file name is expected argument after options\n");
        exit(EXIT_FAILURE);
    }
    inputFile = argv[optind];
    //fprintf(stderr, "The input file name is %s.\n", inputFile);
    
    // Process one file to associate reads with tandem repeats
    time_all = 0; time_memory = 0; time_range = 0; time_period = 0;
    time_initialize_input_string = 0;
    time_count_table = 0;
    time_wrap_around_DP = 0; time_chaining = 0;
    query_counter = 0;
    
    
    struct timeval s, e;
    gettimeofday(&s, NULL);
    
    int read_cnt = handle_one_file(inputFile, print_alignment);
    
    gettimeofday(&e, NULL);
    time_all = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
    if(print_computation_time){
        fprintf(stderr, "Computation time\n");

        fprintf(stderr, "%f\tall\n",           time_all);
        fprintf(stderr, "%f\tallocating memory\n", time_memory);
        fprintf(stderr, "%f\tranges\n",         time_range);
        
        fprintf(stderr, "%f\tComputing periods\n", time_period);
        fprintf(stderr, "\t%f\tInitialize the input\n", time_initialize_input_string);
        fprintf(stderr, "\t%f\tcount table generation\n",   time_count_table);
        fprintf(stderr, "\t%f\twrap around\n",   time_wrap_around_DP);
        fprintf(stderr, "\t%f\tchaining\n",   time_chaining);
        fprintf(stderr, "\t%i\tCount of queries\n", query_counter);
    }
	MPI_Finalize();
    return EXIT_SUCCESS;
}
