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

void search_De_Bruijn_graph( int* a_rep_unit_string, int* rep_unit_score, int query_start, int query_end, int inputLen, int k,
    int *return_predicted_rep_period, int *return_actual_rep_period ){

    // Starting from the initial k-mer, traverse the De Bruijn graph of all k-mers in a greedy manner
    struct timeval s, e;
    gettimeofday(&s, NULL);
    
    // Generate the list of k-mers and search for the node with the maximum count
    init_inputString(k, query_start, query_end, inputLen);
    int width = query_end - query_start + 1;
    int maxNode = generate_freqNode_return_maxNode(query_start, query_end, k, width);
    
    // de Bruijn graph search
    int initialNode = maxNode;
    int Node = initialNode;
    int actual_rep_period = 0;
    for(int i = 0; i < MAX_PERIOD; i++){ a_rep_unit_string[i] = -1; }
    int list_tiebreaks[MAX_tiebreaks];
    int list_tiebreaks_new[MAX_tiebreaks];
    
    for(int l=0; l< MAX_PERIOD; l++){
        a_rep_unit_string[l] = Node / pow4[k-1];  // Memorize the first base
        
        rep_unit_score[l] = freq_node(Node, k, width);    //rep_unit_score[l] = count[Node];
        
        int m, lsd, max_lsd, max_count_lsd, init_tiebreaks, tiebreaks;
        list_tiebreaks[0] = 0;
        init_tiebreaks = 1;
        for(m = 1; m <= k; m++){    // look m bases ahead
            max_lsd = 0;
            max_count_lsd = -1;
            tiebreaks = 0;
            for(int i=0; i < init_tiebreaks; i++){
                for(int j=0; j<4; j++){
                    lsd = 4 * list_tiebreaks[i] + j; // lsd: least significant (4) digit
                    int tmp_count = freq_node( pow4[m] * (Node % pow4[k-m]) + lsd, k, width );
                    // int tmp_count = count[ pow4[m] * (Node % pow4[k-m]) + lsd ];
                    
                    if( max_count_lsd < tmp_count){
                        max_count_lsd = tmp_count;
                        max_lsd = lsd;
                        tiebreaks = 0;
                        list_tiebreaks_new[tiebreaks++] = lsd;
                    }else if(max_count_lsd == tmp_count){ // Allow more tiebreaks
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
        if(Node == initialNode){     // Break if you hit the initial node.
            actual_rep_period = l+1;
            break;
        }
    }
    gettimeofday(&e, NULL);
    time_search_De_Bruijn_graph
    += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
    *return_predicted_rep_period = actual_rep_period;
    *return_actual_rep_period   = actual_rep_period;
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

void polish_repeat(repeat_in_read *rr, int inputLen){
    
    int k = rr->Kmer; // 6 < tt->Kmer
    
    if( rr->rep_period <= k){ // No need to polish the repeat unit
        return;
    }
    init_inputString(k, rr->rep_start, rr->rep_end, inputLen);
    int width = rr->rep_end - rr->rep_start + 1;
    int maxNode = generate_freqNode_return_maxNode(rr->rep_start, rr->rep_end, k, width);
    /*
    memset(count, 0, pow4[k]*4);    //for(int i = 0; i < pow4[k]; i++){ count[i] = 0;}
    for(int i = rr->rep_start; i <= rr->rep_end; i++){
        count[ inputString[i] ]++; }    // Perform counting
    */
    // Convert a char array into an int array
    int int_unit[MAX_PERIOD];
    for(int i = 0; i < rr->rep_period; i++){
        switch (rr->string[i]) {
            case 'A': int_unit[i]=0; break;
            case 'C': int_unit[i]=1; break;
            case 'G': int_unit[i]=2; break;
            case 'T': int_unit[i]=3; break;
            default: fprintf(stderr, "fatal input char %c", rr->string[i]); exit(EXIT_FAILURE);
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
        //int tmp_best_freq = count[refNode];
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


// This function is not effective in increasing the accuracy of predicting repeat units perfectly
void revise_by_progressive_multiple_alignment(
    int* a_rep_unit_string, int rep_period, int query_start, int query_end, int k)
    {
        
    #ifdef DEBUG_progressive_multiple_alignment
        printf("rep_period, query_start, query_end, k = %i %i %i %i\n", rep_period, query_start, query_end, k);
    #endif
        
    int consensus[MAX_PERIOD][5];   // 0-3 for A,C,G,T and 4 for a gap
    int gaps[MAX_PERIOD][4];        // 0-3 for A,C,G,T
        
    for(int i=0; i<MAX_PERIOD; i++){
        for(int j=0; j<5; j++){
            consensus[i][j] = 0;
        }
        for(int j=0; j<4; j++){
            gaps[i][j] = 0;
        }
    }
    
    // A progressive multiple alignment between

    int *matDP = WrapDP;  // Reuse WrapDP by renaming WrapDP by oneDP
    int next = rep_period + 1;
    
    for(int l = query_start; l < query_end; ){
        // Perform a global alignment
        // i scans the reference starting from a_rep_unit_string.
        // j scans from position l in oneInputString.
        
        for(int i=0; i < 2*next*next; i++){ matDP[i] = 0;}
        
        // Initialize the boundaries
        for(int i = 0; i <= rep_period; i++){
            //matDP[ i * next + 0] = 0;
            matDP[ i * next + 0] = (-1) * INDEL_PENALTY * i;
        }
        for(int j = 0; j <= 2*rep_period; j++){
            matDP[ 0 * next + j] = 0;
            //matDP[ 0 * next + j] = (-1) * INDEL_PENALTY * j;
        }
        
        // Inductive step
        for(int i = 0; i < rep_period; i++){
            for(int j = 0; j < rep_period*2 && l+j < query_end; j++ ){
                // Search 2*rep_period positions
                int gain_or_penalty;
                if( a_rep_unit_string[ i ] == orgInputString[ l+j ]){
                    gain_or_penalty = MATCH_GAIN;
                }else{
                    gain_or_penalty = (-1) * MISMATCH_PENALTY;
                }
                     matDP[ (i+1) * next +(j+1)] =
                MAX( matDP[  i    * next + j ]     + gain_or_penalty,
                MAX( matDP[ (i+1) * next + j ]     - INDEL_PENALTY,
                     matDP[  i    * next +(j+1) ]  - INDEL_PENALTY ) );
            }
        }
        
        int max_matDP = 0;
        int max_i = 0;
        int max_j = 0;
        for(int i = 0; i <= rep_period; i++){
            for(int j = 0; j <= 2*rep_period && l+j < query_end; j++ ){
                if( max_matDP < matDP[ i * next + j ] )
                {
                    max_matDP = matDP[ i * next + j ];
                    max_i = i;
                    max_j = j;
                }
            }
        }


#ifdef DEBUG_progressive_multiple_alignment
        printf("l = %i,\tmax_j = %i\n", l, max_j);
#endif
        if(max_j == 0 || max_matDP < 0){ break; }
        // No meaningful alignments are obtained.

        
#ifdef DEBUG_progressive_multiple_alignment
        printf("l, max_j, matDP = %i %i %i\n", l, max_j, max_matDP );
#endif
        
        // traceback and revise the consensus alignment
        int max_val = matDP[ max_i * next + max_j ];
        int i = max_i - 1;
        int j = max_j - 1;
        // max_val is assumed to be matDP[ (i+1) * next + (j+1) ]
        while( i >= 0 && j >= 0){
            if( i >= 0 && j >= 0 &&
                matDP[ i * next + j ] + MATCH_GAIN == max_val)
            {
                consensus[i][orgInputString[ l+j ]]++;
                max_val -= MATCH_GAIN;
                i--; j--;
            }else if( i >= 0 && j >= 0 &&
                matDP[ i * next + j ] - MISMATCH_PENALTY == max_val)
            {
                consensus[i][orgInputString[ l+j ]]++;
                max_val += MISMATCH_PENALTY;
                i--; j--;
            }else if( j >= 0 &&
                matDP[ (i+1) * next + j ] - INDEL_PENALTY  == max_val)
            {   // Insert a base after the i-th position
                gaps[i+1][orgInputString[ l+j ]]++;
                max_val += INDEL_PENALTY;
                j--;
            }else if( i >= 0 &&
                matDP[ i * next + (j+1) ] - INDEL_PENALTY == max_val)
            {   // Delete the base at the i-th position
                consensus[i][4]++;
                max_val += INDEL_PENALTY;
                i--;
            }else{
                fprintf(stderr, "fatal error in progressive multiple alignment DP at (i, j) = (%i %i)\n", i, j);
                exit(EXIT_FAILURE);
            }
        }
        l = l + max_j;
    }
    
    // Obtain the consensus string
    // The current implementation is naive and needs improvement.
    
    int sum_max = 0;
    for(int i=0; i < rep_period; i++){
        int tmp_max = consensus[i][0];
        for(int j=1; j<4; j++){
            if(tmp_max < consensus[i][j]){
                tmp_max = consensus[i][j];
            }
        }
        sum_max += tmp_max;
    }
    double avg_max = (double)sum_max / rep_period;
                
    int revised_rep_unit_string[MAX_PERIOD];

    int revised_i = 0;
    for(int i=0; i < rep_period; i++){
        int max_j = 0;
        int tmp_max = consensus[i][0];
        for(int j=1; j<5; j++){
            if(tmp_max < consensus[i][j]){
                max_j = j;
                tmp_max = consensus[i][j];
            }
        }
        revised_rep_unit_string[revised_i++] = max_j;
        
        int max_gap_j = 0;
        tmp_max = gaps[i][0];
        for(int j=1; j<4; j++){
            if(tmp_max < gaps[i][j]){
                max_gap_j = j;
                tmp_max = gaps[i][j];
            }
        }
        if(avg_max <= tmp_max){
            revised_rep_unit_string[revised_i++] = max_gap_j;
        }
    }
        
#ifdef DEBUG_progressive_multiple_alignment
    printf("rep_period, sum_max, avg_max = %i %i %f\n", rep_period, sum_max, avg_max);
    
    printf("Pred.\t");
    for(int i=0; i<rep_period; i++){
        printf("%i", a_rep_unit_string[i]);
    }
    printf("\nRevised\t");
    for(int i=0; i<revised_i; i++){
        printf("%i", revised_rep_unit_string[i]);
    }
        
    printf("\nconsensus\n");
    for(int j=0; j<5; j++){
        printf("%i\t", j);
        for(int i=0; i<rep_period; i++){
            printf("%i", consensus[i][j]);
        }
        printf("\n");
    }
        printf("gaps\n");
    for(int j=0; j<4; j++){
        printf("%i\t", j);
        for(int i=0; i<rep_period; i++){
            printf("%i", gaps[i][j]);
        }
        printf("\n");
    }
#endif

/*
    for(int i=0; i<revised_i; i++){
        a_rep_unit_string[i] = revised_rep_unit_string[i];
    }
*/
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
            default: fprintf(stderr, "fatal input char %c", string[i]); exit(EXIT_FAILURE);
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

#ifdef DEBUG_progressive_multiple_alignment
    revise_by_progressive_multiple_alignment(int_unit, rep_period, rep_start, rep_end, k);
#endif
}
