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

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

void free_global_variables_and_exit(){
    // If any of the above global variables failed to be allocated in the heap, free other variables and exit.
    if(pow4              != NULL){ free(pow4); }
    if(orgInputString    != NULL){ free(orgInputString); }
    if(inputString       != NULL){ free(inputString); }
    if(inputString_w_rand!= NULL){ free(inputString_w_rand); }
    if(count             != NULL){ free(count); }
    if(sortedString      != NULL){ free(sortedString); }
    if(RRs               != NULL){ free(RRs);}
    if(directional_index_tmp != NULL){ free(directional_index_tmp); }
    if(directional_index != NULL){ free(directional_index); }
    if(directional_index_end != NULL){ free(directional_index_end); }
    if(directional_index_w != NULL){ free(directional_index_w); }
    if(vector0           != NULL){ free(vector0); }
    if(vector1           != NULL){ free(vector1); }
    if(vector2           != NULL){ free(vector2); }
    if(freq_interval_len != NULL){ free(freq_interval_len); }
    //if(count_period_all  != NULL){ free(count_period_all); }
    //if(rep_unit_string   != NULL){ free(rep_unit_string); }
    if(WrapDP            != NULL){ free(WrapDP); }
    if(alignment_input   != NULL){ free(alignment_input); }
    if(alignment_symbols != NULL){ free(alignment_symbols); }
    if(alignment_repeats != NULL){ free(alignment_repeats); }
    for(int i=0; i<(MAX_PERIOD + 1); i++){
        if( consensus[i] != NULL ){ free(consensus[i]); }
    }
    if(consensus         != NULL){ free(consensus); }
    for(int i=0; i<(MAX_PERIOD + 1); i++){
        if( gaps[i] != NULL ){ free(gaps[i]); }
    }
    if( gaps             != NULL ){ free(gaps); }
    fprintf(stderr, "cannot allocate space for one of global variables in the heap.\n");
    exit(EXIT_FAILURE);
}

void malloc_global_variables(){
    
    // Allocate the main memory for global variables in the heap
    pow4  = (int *)malloc(sizeof(int) * (maxKmer+1));
    if( pow4 == NULL ){ free_global_variables_and_exit(); }
    pow4[0] = 1;
    for(int i=1; i<= maxKmer; i++){
        pow4[i] = pow4[i-1] * 4;
    }
    
    orgInputString  = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( orgInputString == NULL ){ free_global_variables_and_exit(); }
    
    inputString     = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( inputString == NULL ){ free_global_variables_and_exit(); }
    
    inputString_w_rand = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( inputString_w_rand == NULL ){ free_global_variables_and_exit(); }
    
    count           = (int *)malloc( sizeof(int) * pow4[maxKmer]);
    if( count == NULL ){ free_global_variables_and_exit(); }
    
    sortedString    = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( sortedString == NULL ){ free_global_variables_and_exit(); }
    
    RRs = (repeat_in_read*) malloc(2*sizeof(repeat_in_read));
    if( RRs == NULL ){ free_global_variables_and_exit(); }
    
    directional_index_tmp   = (double *)malloc(sizeof(double) * MAX_INPUT_LENGTH);
    if( directional_index_tmp == NULL ){ free_global_variables_and_exit(); }
    
    directional_index       = (double *)malloc(sizeof(double) * MAX_INPUT_LENGTH);
    if( directional_index == NULL ){ free_global_variables_and_exit(); }
    
    directional_index_end   = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( directional_index_end == NULL ){ free_global_variables_and_exit(); }
    
    directional_index_w     = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( directional_index_w == NULL ){ free_global_variables_and_exit(); }
    
    vector0 = (int *)malloc(sizeof(int) * 4 * BLK);
    if( vector0 == NULL ){ free_global_variables_and_exit(); }
    
    vector1 = (int *)malloc(sizeof(int) * 4 * BLK);
    if( vector1 == NULL ){ free_global_variables_and_exit(); }
    
    vector2 = (int *)malloc(sizeof(int) * 4 * BLK);
    if( vector2 == NULL ){ free_global_variables_and_exit(); }
    
    freq_interval_len   = (double *)malloc( sizeof(double) * MAX_INPUT_LENGTH);
    if( freq_interval_len == NULL ){ free_global_variables_and_exit(); }
    
    //count_period_all    = (int *)malloc( sizeof(int) * MAX_PERIOD);
    //if( count_period_all == NULL ){ free_global_variables_and_exit(); }
    
    //rep_unit_string     = (int *)malloc( sizeof(int) * MAX_PERIOD);
    //if( rep_unit_string == NULL ){ free_global_variables_and_exit(); }
    
    WrapDP          = (int *)malloc(sizeof(int) * WrapDPsize);
    if( WrapDP == NULL ){ free_global_variables_and_exit(); }
    
    alignment_input = (char *)malloc(sizeof(char) * MAX_INPUT_LENGTH);
    if( alignment_input == NULL ){ free_global_variables_and_exit(); }
    
    alignment_symbols = (char *)malloc(sizeof(char) * MAX_INPUT_LENGTH);
    if( alignment_symbols == NULL ){ free_global_variables_and_exit(); }
    
    alignment_repeats = (char *)malloc(sizeof(char) * MAX_INPUT_LENGTH);
    if( alignment_repeats == NULL ){ free_global_variables_and_exit(); }
    
    consensus = malloc(sizeof(int *) * (MAX_PERIOD + 1));
    if( consensus == NULL ){ free_global_variables_and_exit(); }
    for(int i=0; i<(MAX_PERIOD + 1); i++){
        consensus[i] = malloc(sizeof(int) * 5);
        if( consensus[i] == NULL ){ free_global_variables_and_exit(); }
    }
    
    gaps = malloc(sizeof(int *) * (MAX_PERIOD + 1));
    if( gaps == NULL ){ free_global_variables_and_exit(); }
    for(int i=0; i<(MAX_PERIOD + 1); i++){
        gaps[i] = malloc(sizeof(int) * 4);
        if( gaps[i] == NULL ){ free_global_variables_and_exit(); }
    }
}

void free_global_variables(){
    free(pow4);
    free(orgInputString);
    free(inputString);
    free(inputString_w_rand);
    free(count);
    free(sortedString);
    free(RRs);
    free(directional_index_tmp);
    free(directional_index);
    free(directional_index_end);
    free(directional_index_w);
    free(vector0);
    free(vector1);
    free(vector2);
    free(freq_interval_len);
    //free(count_period_all);
    //free(rep_unit_string);
    free(WrapDP);
    free(alignment_input);
    free(alignment_symbols);
    free(alignment_repeats);
    for(int i=0; i<(MAX_PERIOD + 1); i++){ free(consensus[i]); }
    for(int i=0; i<(MAX_PERIOD + 1); i++){ free(gaps[i]); }
    free(consensus);
    free(gaps);
    
}

int handle_one_file(char *inputFile, int print_multiple_TR, int print_alignment){
    //---------------------------------------------------------------------------
    // Feed a string from a file, convert the string into a series of integers
    //---------------------------------------------------------------------------
    
    FILE *fp = fopen(inputFile, "r");
    if(fp == NULL){
        fprintf(stderr, "fatal error: cannot open %s\n", inputFile);
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
    
    char s[BLK];
    char readID[BLK];
    int i, inputLen, charCode;
    int cnt=0;
    int read_cnt = 0;
    int firstRead = 1;  // 1 means the first read.
    
    struct timeval s_time, e_time;
    gettimeofday(&s_time, NULL);
    malloc_global_variables();
    gettimeofday(&e_time, NULL);
    time_memory += (e_time.tv_sec - s_time.tv_sec) + (e_time.tv_usec - s_time.tv_usec)*1.0E-6;

    
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
                handle_one_read(readID, inputLen, read_cnt, print_multiple_TR, print_alignment);
                read_cnt++;
            }
            // Feed the header of the read.
            for(i=1; s[i]!='\0' && s[i]!='\n' && i<BLK; i++){
                readID[i-1] = s[i];
            }
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
    handle_one_read(readID, inputLen, read_cnt, print_multiple_TR, print_alignment);
    read_cnt++;

    
    fclose(fp);

    free_global_variables();
    
    return(read_cnt);
}
