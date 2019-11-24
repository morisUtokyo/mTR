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

void init_inputString(int k, int query_start, int query_end, int inputLen){
    
    struct timeval s, e;
    gettimeofday(&s, NULL);
    
    for(int i = query_start; i < query_end+k-1 && i < inputLen; i++){
        inputString[i] = orgInputString[i];
    }
    int tmp = 0;

    for(int i = query_start; i < query_start+k-1; i++){
        tmp = 4 * tmp + inputString[i];  // compute 4 decimal of first k-1 letters
    }
    for(int i = query_start; i < query_end && i < (inputLen-k+1); i++){
        inputString[i] = 4 * tmp + inputString[i+k-1];
        tmp = inputString[i] % pow4[k-1]; //　remainder, compute 4 decimal of length k-1
        if(tmp < 0){
            fprintf(stderr, "fatal error at %i\t %i \t%i\n", i, inputString[i], pow4[k-1] );
            exit(EXIT_FAILURE);
        }
    }
    gettimeofday(&e, NULL);
    time_initialize_input_string += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
}

int nearPrime(int width){
    
    int Primes[9]={1009, 2003, 4001, 8009, 16001, 32003, 64007, 128021, 256019};
    int i=0;
    for(i=0; Primes[i] < width && i<8; i++ );
    return(Primes[i+1]);    // if width <= 2003, return 4001
    
}

int generate_freqNode_return_maxNode(int query_start, int query_end, int k, int width){
    
    struct timeval s, e;
    gettimeofday(&s, NULL);
    
    int Prime = nearPrime(width);
    
    int maxNode;
    
    if(k <= count_maxKmer){     // 4^k
        memset(count, 0, pow4[k]*4);    //for(int i = 0; i < pow4[k]; i++){ count[i] = 0;}
        for(int i = query_start; i <= query_end; i++){ count[ inputString[i] ]++; }
        int max_count = -1;
        for(int i = query_start; i <= query_end; i++){
            if( max_count < count[ inputString[i] ] ){
                max_count = count[ inputString[i] ];
                maxNode = inputString[i];
            }
        }
    }else{
        for(int i = 0; i < Prime; i++){
            freqNode[i][0] = -1;    // no entry of Node
            freqNode[i][1] = 0;     // reset the frequency to 0
        }

        int maxFreq = -1;
        int Node;
        for(int i = query_start; i <= query_end; i++){
            Node = inputString[i];
            int h =  Node % Prime;     // a hash value of the node
            for(int j=0; ; j++){
                // j counts the number of trials
                if(freqNode[h][0] == -1){
                    // Found an empty entry
                    freqNode[h][0] = Node;
                    freqNode[h][1] = 1;
                    break;
                }else if(freqNode[h][0] == Node){
                    // The node has been already put into the table. Increment the frequency.
                    freqNode[h][1]++;
                    break;
                }else if(Prime <= j){
                    fprintf(stderr, "The hash table for nodes is full.\n");
                    exit(EXIT_FAILURE);
                }else{
                    // Move to the next entry
                    h = (h+1) % Prime;
                }
            }
            if(maxFreq < freqNode[h][1]){
                maxFreq = freqNode[h][1];
                maxNode = Node;
            }
        }
    }
    gettimeofday(&e, NULL);
    time_count_table += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
    return(maxNode);
}

int freq_node(int Node, int k, int width){
    
    if(k <= count_maxKmer){
        return( count[Node] );
    }else{
        int Prime = nearPrime(width);
        int h = Node % Prime;
        for(int j=0;  ; j++){
            if(freqNode[h][0] == Node){
                return(freqNode[h][1]);
            }else if(freqNode[h][0] == -1){
                // The node is absent in the table.
                return(0);
            }else if(Prime <= j){
                fprintf(stderr, "The hash table for nodes is full.\n");
                exit(EXIT_FAILURE);
            }else{
                // Move to the next
                h = (h+1) % Prime;
            }
        }
    }
}


void print_4_decimal_array(int* val, int len, char *return_string){
    strcpy(return_string, "");
    for(int i=0; i<len; i++){
        switch(val[i]){
            case 0: strcat(return_string, "A"); break;
            case 1: strcat(return_string, "C"); break;
            case 2: strcat(return_string, "G"); break;
            case 3: strcat(return_string, "T"); break;
            default: fprintf(stderr, "fatal error: input char %i", val[i]); exit(EXIT_FAILURE);
        }
    }
}

void search_De_Bruijn_graph( int query_start, int query_end, repeat_in_read *rr){
    // Starting from the initial k-mer, traverse the De Bruijn graph of all k-mers in a greedy manner
    int rep_unit_string[MAX_PERIOD];
    int rep_unit_score[MAX_PERIOD];
    int inputLen = rr->inputLen;
    int k = rr->Kmer;
    
    // Generate the list of k-mers and search for the node with the maximum count
    struct timeval s, e;
    gettimeofday(&s, NULL);
    init_inputString(k, query_start, query_end, inputLen);
    int width = query_end - query_start + 1;
    int maxNode = generate_freqNode_return_maxNode(query_start, query_end, k, width);
    gettimeofday(&e, NULL);
    time_init_search_De_Bruijn_graph += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
    // de Bruijn graph search
    int initialNode = maxNode;
    int Node = initialNode;
    int actual_rep_period = 0;
    for(int i = 0; i < MAX_PERIOD; i++){ rep_unit_string[i] = -1; }
    int list_tiebreaks[MAX_tiebreaks];
    int list_tiebreaks_new[MAX_tiebreaks];
    
    for(int l=0; l < MAX_PERIOD && l < (query_end - query_start)/MIN_NUM_FREQ_UNIT; l++){
        // l must be smaller than the range of the query
        rep_unit_string[l] = Node / pow4[k-1];  // Memorize the first base
        rep_unit_score[l] = freq_node(Node, k, width);    //rep_unit_score[l] = count[Node];
        
        int m, lsd, max_lsd, max_count_lsd, init_tiebreaks, tiebreaks, max_lookahead;
        list_tiebreaks[0] = 0;
        init_tiebreaks = 1;
        // Look k bases ahead if the period is of >10 bases.
        if( l < 10 ){
            max_lookahead = 1;
        }else{
            max_lookahead = k;
        }
        for(m = 1; m <= max_lookahead; m++){    // look m bases ahead
            max_lsd = 0;
            max_count_lsd = -1;
            tiebreaks = 0;
            for(int i=0; i < init_tiebreaks; i++){
                for(int j=0; j<4; j++){
                    lsd = 4 * list_tiebreaks[i] + j; // lsd: least significant (4) digit
                    int tmp_count = freq_node( pow4[m] * (Node % pow4[k-m]) + lsd, k, width );
                    if( max_count_lsd < tmp_count){
                        max_count_lsd = tmp_count;
                        max_lsd = lsd;
                        tiebreaks = 0;
                        list_tiebreaks_new[tiebreaks++] = lsd;
                    }else if(max_count_lsd == tmp_count){ // Allow more tiebreaks
                        //max_lsd = lsd;
                        if(tiebreaks < MAX_tiebreaks){
                            list_tiebreaks_new[tiebreaks++] = lsd;
                        }
                    }
                }
            }
            if( tiebreaks == 1 ){    // Break if no tiebreak
                break;
            }else{
                for(int i=0; i<tiebreaks; i++){
                    list_tiebreaks[i] = list_tiebreaks_new[i];
                }
                init_tiebreaks = tiebreaks;
            }
        }
        Node = 4 * (Node % pow4[k-1]) + (max_lsd / pow4[m-1] ); // Shift by one base
        if(Node == initialNode){
            actual_rep_period = l+1;
            break;
        }
    }
    if(actual_rep_period < MIN_PERIOD){ // if actual_rep_period == 0, return NULL.
        return;
    }else{
        rr->rep_period           = actual_rep_period;
        rr->predicted_rep_period = actual_rep_period;
        rr->ConsensusMethod      = DeBruijnGraphSearch;
        print_4_decimal_array(rep_unit_string, actual_rep_period, rr->string);
        for(int i=0; i<rr->rep_period; i++){
            rr->string_score[i] = rep_unit_score[i];
        }
        freq_2mer_array(rep_unit_string, actual_rep_period, rr->freq_2mer);
    }

}

int score_for_alignment(int start, int k, int bestNode, int rep_period, int* int_unit, int width){
    
    int tmpNode = bestNode;
    int sumFreq = 0;
    for( int j = start; start - k < j; j--){ // wrap around
        tmpNode = int_unit[ j%rep_period ] * pow4[k-1] + (tmpNode / 4) ;
        sumFreq += freq_node(tmpNode, k, width);
        //sumFreq += count[tmpNode];
    }
    return(sumFreq);
}

int suspicious(repeat_in_read *rr, int j){
    // More than 80% of the rr->Kmer bases before j are covered by only one base.
    int cnt = 0;
    for(int i=0; i < (rr->Kmer - 1) && 0 <= j-i; i++){
        if(rr->string_score[j-i] < 2){ cnt++; }
    }
    if( (rr->Kmer-1)*0.8 < (double)cnt ){
        return(1);
    }else{
        return(0);
    }
}

void polish_repeat(repeat_in_read *rr){
    
    int k = rr->Kmer;
    int inputLen = rr->inputLen;
    
    if( rr->rep_period <= k){ // No need to polish the repeat unit
        return;
    }
    init_inputString(k, rr->rep_start, rr->rep_end, inputLen);
    int width = rr->rep_end - rr->rep_start + 1;
    int maxNode = generate_freqNode_return_maxNode(rr->rep_start, rr->rep_end, k, width);

    // Convert a char array into an int array
    int int_unit[MAX_PERIOD];
    for(int i = 0; i < rr->rep_period; i++){
        switch (rr->string[i]) {
            case 'A': int_unit[i]=0; break;
            case 'C': int_unit[i]=1; break;
            case 'G': int_unit[i]=2; break;
            case 'T': int_unit[i]=3; break;
            default: fprintf(stderr, "polish_repeat: fatal input char %c\n", rr->string[i]); exit(EXIT_FAILURE);
        }
    }
    int rep_period = rr->rep_period;
    int int_revised_unit[MAX_PERIOD];
    int j_revised = MAX_PERIOD - 1;
    
    int refNode = 0;
    for(int i = 0; i < k; i++){
        refNode = int_unit[i] * pow4[k-1-i] + refNode;
        // Treat the first k bases as the end of the input string
    }
    
    int bestNode = refNode;
    for(int j = rep_period-1; 0 <= j; ){
        int refNode = int_unit[j] * pow4[k-1] + (bestNode / 4);
        int tmp_best_freq = freq_node(refNode, k, width);
        bestNode = refNode;
        
        if(rr->string_score[j] == 1 && suspicious(rr, j) == 1 )
        {
            for(int l = 0; l < 4; l++){
                int alternativeNode = refNode + (l - int_unit[j]) * pow4[k-1];
                if( alternativeNode < 0 || pow4[k] <= alternativeNode){
                    fprintf(stderr, "fatal error: alternativeNode is out of range.\n");
                    exit(EXIT_FAILURE);
                }
                    // Replace the first base with l = 0, 1, 2, 3.
                if( tmp_best_freq < freq_node(alternativeNode, k, width) ){
                    tmp_best_freq = freq_node(alternativeNode, k, width);
                    bestNode = alternativeNode;
                }
            }
            if(bestNode == refNode){    // No need to revise the next node
                int_revised_unit[ j_revised-- ] = int_unit[ j-- ];
            }else{  // Use next three bases to calculate the score
                int score_del = score_for_alignment( j, k, bestNode, rep_period, int_unit, width);
                int score_sub = score_for_alignment( j-1, k, bestNode, rep_period, int_unit, width);
                int score_ins = -1;
                if( bestNode / pow4[k-1] == int_unit[ (j-1) % rep_period] ){
                    score_ins = score_for_alignment( j-2, k, bestNode, rep_period, int_unit, width);
                }
                int_revised_unit[ j_revised-- ] = bestNode / pow4[k-1]; // The first base

                int max_score = MAX( MAX(score_del, score_sub), score_ins);
                if( max_score == score_del){
                    j = j;  // Do not increment j as int_unit[j] is used in the next step.
                }else if(max_score == score_sub){
                    j -= 1;
                }else{
                    j = j-2;
                }
            }
        }else{
            int_revised_unit[ j_revised-- ] = int_unit[ j-- ];
        }
        if(j_revised < 0){ // fails to revise
            return;
        }
    }
    rr->rep_period = (MAX_PERIOD-1) - j_revised;
    // Change Num_matches etc. using wrap around DP
    for(int i = 0; i < rr->rep_period; i++){
        char c;
        switch( int_revised_unit[ i + j_revised + 1 ] ) {
            case 0: c = 'A'; break;
            case 1: c = 'C'; break;
            case 2: c = 'G'; break;
            case 3: c = 'T'; break;
            default: fprintf(stderr, "fatal input int %i at %i\n", int_revised_unit[i+j_revised+1], i ); exit(EXIT_FAILURE);
        }
        rr->string[i] = c;
    }
    rr->string[rr->rep_period] = '\0';
}

// Columns from the top represent mismatch ratios: 0.25, 0.2, 0.15, 0.1, 0.05, and 0.02.
// Rows are read coverages that range from 0 to 20.
int min_missing_bases[6][21] = {
    0,1,2,2,3,3,3,4,4,4,5,5,6,6,6,7,7,7,8,8,8,
    0,1,1,2,2,3,3,3,4,4,4,5,5,5,5,6,6,6,7,7,7,
    0,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,
    0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,
    0,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,
    0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2
};

int min_missing(double ratio, int coverage){
    int col, row;
    
    if(ratio > 0.225){        row = 0; }
    else if(ratio > 0.175 ){  row = 1; }
    else if(ratio > 0.125 ){  row = 2; }
    else if(ratio > 0.075 ){  row = 3; }
    else if(ratio > 0.025 ){  row = 4; }
    else{                     row = 5; }
    
    col = coverage;
    
    return(min_missing_bases[row][col]);
}

// This function is not effective in increasing the accuracy of predicting repeat units perfectly
//void revise_representative_unit( repeat_in_read *rr,  char *string, int unit_len, int query_start, int query_end){
void revise_representative_unit( repeat_in_read *rr){

    // Initialization
    char *string    = rr->string;
    int unit_len    = rr->rep_period;
    int query_start = rr->rep_start;
    int query_end   = rr->rep_end;
    
    // Convert a char array into an int array
    int rep_unit[MAX_PERIOD];
    for(int i = 0; i < unit_len; i++){
        switch (string[i]) {    // Shift by 1
            case 'A': rep_unit[i+1] = 0; break;
            case 'C': rep_unit[i+1] = 1; break;
            case 'G': rep_unit[i+1] = 2; break;
            case 'T': rep_unit[i+1] = 3; break;
            default: fprintf(stderr, "revise rep: fatal input char %c\n", string[i]); exit(EXIT_FAILURE);
        }
    }
    
    int *rep;
    rep = &orgInputString[query_start];
    int rep_len = query_end - query_start + 1;
    
    int i, j;
    int next = unit_len+1;

    for(j=0; j<=rep_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = 0; // local alignment
    }
    
    int max_wrd = 0;
    int max_i = 0;
    int max_j = 0;
    int val_match, val_mismatch, val_insertion, val_deletion;
    for(i=1; i <= rep_len; i++){      // Scan repeat
        for(j=1; j<=unit_len; j++){   // Scan rep_unit
            if( WrapDPsize <= next*i + j ){
                fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
                exit(EXIT_FAILURE);
            }
            if(rep[i] == rep_unit[j]){    // *1*-origin index !!!!
                WrapDP[next*i + j] = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
            }else{
                val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
                val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
                if(j > 1){
                    val_deletion = WrapDP[next*i + j-1] - INDEL_PENALTY;
                    WrapDP[next*i + j] = MAX(0, MAX( MAX( val_mismatch, val_insertion), val_deletion));
                }else{
                    WrapDP[next*i + j] = MAX(0, MAX( val_mismatch, val_insertion));
                }
            }
            if(max_wrd < WrapDP[next*i + j])
            {
                max_wrd = WrapDP[next*i + j];
                max_i = i;
                max_j = j;
            }
        }
        // wrap around
        WrapDP[next*i + 0] = WrapDP[next*i + unit_len];
    }
    
    // trace back the optimal alignment while storing it in the data structure "alignment"
    int consensus[MAX_PERIOD][5];   // 0-3 for A,C,G,T and 4 for a gap
    int missing[MAX_PERIOD][4];        // 0-3 for A,C,G,T
        
    for(int j=0; j<MAX_PERIOD; j++){
        for(int q=0; q<5; q++){ consensus[j][q] = 0;   }
        for(int q=0; q<4; q++){ missing[j][q] = 0;    }
    }
    i = max_i;
    j = max_j;
    if(j == 0){ j = unit_len; } // 1-origin index
    // global alignment
    while(i > 0 && WrapDP[next*i + j] > 0){
        val_match       = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
        val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
        val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
        val_deletion    = WrapDP[next*i + j-1]      - INDEL_PENALTY;
        
        if( max_wrd == val_match          && rep[i] == rep_unit[j]){
            consensus[j][rep[i]]++; // Increment the number of match rep[i]
            max_wrd -= MATCH_GAIN;
            i--; j--;
        }else if( max_wrd == val_mismatch && rep[i] != rep_unit[j]){     // mismatch
            consensus[j][rep[i]]++; // Increment the number of mismatch rep[i]
            max_wrd += MISMATCH_PENALTY;
            i--; j--;
        }else if( max_wrd == val_deletion){
            consensus[j][4]++;      // Increment the number of deletions from the representative
            max_wrd += INDEL_PENALTY;
            j--;
        }else if( max_wrd == val_insertion){    // insertion
            missing[j][rep[i]]++;      // Increment the number of the base missing from the representative
            max_wrd += INDEL_PENALTY;
            i--;
        }else if( max_wrd == 0){
            break;
        }else{
            fprintf(stderr, "fatal error in wrap-around DP max_wrd = %i\n", max_wrd);
            exit(EXIT_FAILURE);
        }
        if(j == 0){
            j = unit_len;
        }
    }
    
    int revised_rep_unit[MAX_PERIOD];
    int revised_rep_j = 0;  // 0-origin index
    
    int rep_unit_before[MAX_PERIOD];
    int rep_unit_after[MAX_PERIOD];
    int rep_j = 1;          // 1-origin index
    for(int j=1; j<=unit_len; j++){ // 1-origin index
        int max_v = -1;
        int max_base = -1;
        for(int q=0; q<5; q++){
            if(max_v < consensus[j][q]){
                max_v = consensus[j][q];
                max_base = q;
            }
        }
        rep_unit_before[rep_j] = rep_unit[j];
        rep_unit_after[rep_j] = max_base;
        rep_j++;
        
        if(max_base < 4){ // Non-gap
            revised_rep_unit[revised_rep_j++] = max_base;
        }

        max_v = -1;
        int max_missing  = -1;
        for(int q=0; q<4; q++){  // Assume that missing bases are between the first and last bases
            if(max_v < missing[j][q]){  // If the max value is greater than the previous base.
                max_v = missing[j][q];
                max_missing = q;
            }
        }
        //********************************************************
        int coverage = rr->repeat_len / rr->rep_period;
        if( 5 <= coverage && coverage <= 20 ){
            double mismatch_ratio = (double)(rr->Num_mismatches + rr->Num_insertions + rr->Num_deletions) / rr->repeat_len;
            //double indel_ratio = (double)(rr->Num_insertions + rr->Num_deletions) / rr->repeat_len;
            if( min_missing(mismatch_ratio, coverage) <= max_v ){
                rep_unit_before[rep_j] = 4;
                rep_unit_after[rep_j] = max_missing;
                rep_j++;
                revised_rep_unit[revised_rep_j++] = max_missing;
            }
        }
    }
    print_4_decimal_array(revised_rep_unit, revised_rep_j, rr->string);
    
    
#ifdef DEBUG_revise_representative_unit
    // revise_representative_unit(string, rep_period, rep_start, rep_end);
    printf("\nRevision by multiple alignment\nbef\t");
    for(int i = 1; i < rep_j; i++){
        printf("%i", rep_unit_before[i]);
    }
    printf("\naft\t");
    for(int i = 1; i < rep_j; i++){
        printf("%i", rep_unit_after[i]);
    }
    printf("\nmat\t");
    for(int i = 1; i < rep_j; i++){
        if( rep_unit_before[i] == rep_unit_after[i] ){
            printf(" ");
        }else{
            if(rep_unit_before[i] == 4){
                printf("+");
            }else if(rep_unit_after[i] == 4){
                printf("-");
            }else{
                printf("x");
            }
        }
    }
    printf("\nrev\t");
    for(int i = 0; i < revised_rep_j; i++){     // 0-origin index
        printf("%i", revised_rep_unit[i]);
    }
    printf("\n");
#endif
    
}



void print_freq(int rep_start, int rep_end, int rep_period, char* string, int inputLen, int k){
    
    init_inputString(k, rep_start, rep_end, inputLen);
    int width = rep_end - rep_start + 1;
    int maxNode = generate_freqNode_return_maxNode(rep_start, rep_end, k, width);
    /*
    memset(count, 0, pow4[k]*4); //for(int i = 0; i < pow4[k]; i++){ count[i] = 0;}
    for(int i = rep_start; i <= rep_end; i++){
        count[ inputString[i] ]++; }    // Perform counting
    */
    
    // Convert a char array into an int array
    int int_unit[MAX_PERIOD];
    for(int i = 0; i < rep_period; i++){
        switch (string[i]) {
            case 'A': int_unit[i]=0; break;
            case 'C': int_unit[i]=1; break;
            case 'G': int_unit[i]=2; break;
            case 'T': int_unit[i]=3; break;
            default: fprintf(stderr, "print_freq: fatal input char %c\n", string[i]); exit(EXIT_FAILURE);
        }
    }
    int int_Nodes[MAX_PERIOD];
    int tmp = 0;
    for(int i = 0; i < k-1; i++){
        tmp = 4 * tmp + int_unit[i];    // tmp is a string of length k-1
    }
    for(int i = 0; i < rep_period; i++){
        int_Nodes[i] = 4 * tmp + int_unit[ (i+k-1) % rep_period ];
        tmp = int_Nodes[i] % pow4[k-1]; // remove the first base
    }
    for(int i = 0; i < rep_period; i++){
        int freq = freq_node( int_Nodes[i], k, width);
        //int freq = count[int_Nodes[i]];
        if(freq<10){
            printf("%i", freq);
        }else{
            printf("*");
        }
    }
    printf("\n");

}
