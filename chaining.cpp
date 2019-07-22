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
#include "chaining.h"
#include "mTR.h"
#include <stdio.h>
#include <iostream>
#include <set>
#include <map>
#include <string.h>
using namespace std;
//#define DEBUG_chaining

#define MH_distance_threshold 0.3  // Two 2mer frequency distributions are identical if their Manhattan distance is less than or equal to this threshold.  A small threshold generates smaller groups of repeat units.
#define MIN_REP_LEN_CONCATENATE 100
#define MIN_DISTANCE_CONCATENATE 2000

class Alignment{
public:
	int start_x, start_y, end_x, end_y;
	int initial_score, score;
    Alignment* predecessor;
    
    char    readID[BLK];
    int     inputLen;
    int     rep_start;
    int     rep_end;
    int     repeat_len;
    int     rep_period;
    int     Num_freq_unit;
    int     Num_matches;
    int     Num_mismatches;
    int     Num_insertions;
    int     Num_deletions;
    int     Kmer;
    int     ConsensusMethod;     // 0 = progressive multiple alignment, 1 = De Bruijn graph search
    char    string[MAX_PERIOD];

    Alignment(char* a_readID,
              int   a_inputLen,
              int   a_rep_start,
              int   a_rep_end,
              int   a_repeat_len,
              int   a_rep_period,
              int   a_Num_freq_unit,
              int   a_Num_matches,
              int   a_Num_mismatches,
              int   a_Num_insertions,
              int   a_Num_deletions,
              int   a_Kmer,
              int   a_ConsensusMethod,
              char* a_string )
    {
		start_x     = a_rep_start;
		start_y	    = a_rep_start;
		end_x	    = a_rep_end;
		end_y	    = a_rep_end;
		initial_score   = a_Num_matches;
		score	        = a_Num_matches;
		predecessor	= NULL;
        
        strcpy( readID, a_readID );
        inputLen    = a_inputLen;
        rep_start   = a_rep_start;
        rep_end     = a_rep_end;
        repeat_len  = a_repeat_len;
        rep_period  = a_rep_period;
        Num_freq_unit   = a_Num_freq_unit;
        Num_matches     = a_Num_matches;
        Num_mismatches  = a_Num_mismatches;
        Num_insertions  = a_Num_insertions;
        Num_deletions   = a_Num_deletions;
        Kmer        = a_Kmer;
        ConsensusMethod = a_ConsensusMethod;     // 0 = progressive multiple alignment, 1 = De Bruijn graph search
        strcpy( string, a_string );
        
	}
    
    void print_one_TR(int print_alignment)const{
        printf(
               "%.50s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%s\n",
               //rr.ID,
               readID,
               inputLen,
               rep_start,
               rep_end,
               repeat_len,
               rep_period,
               Num_freq_unit,
               Num_matches,
               (float)Num_matches/repeat_len,
               Num_mismatches,
               Num_insertions,
               Num_deletions,
               Kmer,
               ConsensusMethod,
               string
               );

#ifdef Overlapping
        if(predecessor != NULL){
            if(predecessor->rep_end  >= rep_start){
                printf("------------ overlapping ----------------\n");
            }
        }
#endif
        
        if(print_alignment == 1){
            printf("\n");
            pretty_print_alignment( const_cast<char*>(string), rep_period, rep_start, rep_end);
        }
        
        fflush(stdout);
        
    }
    void print_all_TRs(int print_alignment)const{
        if(predecessor != NULL){
            predecessor->print_all_TRs(print_alignment);
        }
        print_one_TR(print_alignment);
    }
	void print_one_alignment()const{
		cout << "(" <<	start_x << "," << start_y << ")\t-> (" << end_x << ","  << end_y <<
                        ")\t units = " << rep_period << " x " << Num_freq_unit <<
                        "\tscores = " << initial_score << "\t" << score << endl;
	}
    void print_chain()const{
        if(predecessor != NULL){
            predecessor->print_chain();
        }
        print_one_alignment();
    }
	bool isStart(int x_coord)const{
		if(x_coord == start_x)
			return true;
		else
			return false;
	}
	void set_predecessor(Alignment* a){
		predecessor = a;
		score += a->score;
	}
    
    int char2int(char c){
        switch(c){
            case 'A': return(0);
            case 'C': return(1);
            case 'G': return(2);
            case 'T': return(3);
            default: fprintf(stderr, "fatal input char %c", c);
        }
    }
    
    void freq_2mer_array(char* st, int *freq_2mer){
        //void freq_2mer_array(char* st, int len, int *freq_2mer){
        for(int i=0; i<16; i++){
            freq_2mer[i] = 0;
        }
        int len;
        for(len = 0; st[len] != '\0'; len++);
        int index;
        for(int i=1; i<len; i++){
            index = char2int(st[i-1]) * 4 + char2int(st[i]);
            freq_2mer[index]++;
        }
        // wrap around and concatenate the last and first characters
        index = char2int(st[len-1]) * 4 + char2int(st[0]);
        freq_2mer[index]++;
    }
    
    int TRs_in_neighborhood(char* a_string, char* query_string){
        //int TRs_in_neighborhood(int *a_freq_2mer, int *query_freq_2mer){
        // return +1 if they are in neighborhood, and -1 otherwise
        
        int a_freq_2mer[16];
        int query_freq_2mer[16];
        freq_2mer_array(a_string,       a_freq_2mer);
        freq_2mer_array(query_string,   query_freq_2mer);
        
        int diff = 0;
        int query_unit_len = 0;
        for(int i=0; i<16; i++){
            int tmp_diff =  a_freq_2mer[i] - query_freq_2mer[i];
            if(tmp_diff < 0){ tmp_diff = (-1) * tmp_diff; }
            diff += tmp_diff;
            query_unit_len += query_freq_2mer[i];
        }
        // We here set the threshold to a temporary value, but it should be revised.
        query_unit_len++;
        if( diff < MH_distance_threshold * query_unit_len )
            return(1);
        else
            return(-1);
    }
    
    Alignment* extend_alignment(int focal_start, char* focal_string){
        
        if(focal_start - rep_end < MIN_DISTANCE_CONCATENATE){ // Within the distance of 1000bp
            if(TRs_in_neighborhood(string, focal_string) == 1){
                // Strings are similar.
                if(predecessor != NULL){
                    Alignment* answer = predecessor->extend_alignment(rep_start, focal_string);
                    // Start from this alignment.
                    if(answer == NULL){     // No maximum alignment ahead, and answer this.
                        return this;
                    }else{                  // Anoswer the maximum alignment ahread.
                        return answer;
                    }
                }else{
                    return this;
                }
            }else{
                if(predecessor != NULL){
                    return predecessor->extend_alignment(focal_start, focal_string);
                }else{
                    return NULL;
                }
            }
        }else{
            return NULL;
        }
    }
    
    void concatenate_similar_alignments(){
        if(repeat_len > MIN_REP_LEN_CONCATENATE){  // qualified
            if(predecessor != NULL){
                Alignment* maximimal = predecessor->extend_alignment(rep_start, string);
                if(maximimal != NULL)
                {
                    //fprintf(stderr, "Pair of alignments-------------------\n");
                    //print_one_TR(0);
                    //maximimal->print_one_TR(0);
                    
                    int a_rep_unit[MAX_PERIOD];
                    int actual_rep_period = rep_period;
                    int query_start = maximimal->rep_start;
                    int query_end   = rep_end;
                    
                    for(int i=0; i<rep_period; i++){ a_rep_unit[i] = char2int(string[i]); }
                    
                    //init_inputString(Kmer, query_start, query_end, inputLen);
                    //actual_rep_period = search_De_Bruijn_graph(a_rep_unit, query_start, query_end, inputLen, Kmer);
                    
                    wrap_around_DP(a_rep_unit, actual_rep_period,
                                   query_start, query_end,
                                   &rep_start, &rep_end, &repeat_len, &Num_freq_unit,
                                   &Num_matches, &Num_mismatches, &Num_insertions, &Num_deletions);
                    predecessor = maximimal->predecessor;
                    //print_one_TR(0);
                    
                    // Search for another possible concatenation
                    if(maximimal->predecessor != NULL){
                        maximimal->predecessor->concatenate_similar_alignments();
                    }
                }else{
                    predecessor->concatenate_similar_alignments();
                }
            }
        }else{
            if(predecessor != NULL){
                predecessor->concatenate_similar_alignments();
            }
        }
    }
};

set<Alignment*> set_of_alignments;

void insert_an_alignment_into_set(
                      char* readID,
                      int   inputLen,
                      int   rep_start,
                      int   rep_end,
                      int   repeat_len,
                      int   rep_period,
                      int   Num_freq_unit,
                      int   Num_matches,
                      int   Num_mismatches,
                      int   Num_insertions,
                      int   Num_deletions,
                      int   Kmer,
                      int   ConsensusMethod,
                      char* string )
{
    set_of_alignments.insert(
                      new Alignment(
                            readID,
                            inputLen,
                            rep_start,
                            rep_end,
                            repeat_len,
                            rep_period,
                            Num_freq_unit,
                            Num_matches,
                            Num_mismatches,
                            Num_insertions,
                            Num_deletions,
                            Kmer,
                            ConsensusMethod,
                            string ));
}

void search_max(int print_alignment){
    if(set_of_alignments.empty()){
        return;
    }
    int max_matches;
    Alignment* max_alignment;
    // Sort alignments with their start and end positions
    for(set<Alignment*>::iterator iter = set_of_alignments.begin();
        iter != set_of_alignments.end();
        iter++)
    {
        if(iter == set_of_alignments.begin()){
            max_matches = (*iter)->Num_matches;
            max_alignment = (*iter);
        }else if(max_matches < (*iter)->Num_matches){
            max_matches = (*iter)->Num_matches;
            max_alignment = (*iter);
        }
    }
    max_alignment->print_one_TR(print_alignment);
    
    for(set<Alignment*>::iterator iter = set_of_alignments.begin();
        iter != set_of_alignments.end(); iter++){
        set_of_alignments.erase(iter);
    }
}


void chaining(int print_alignment){
    if(set_of_alignments.empty()){
        return;
    }
    multimap<int, Alignment*> sorted_by_X;
    multimap<int, Alignment*> sorted_by_Y;
    
    // Sort alignments with their start and end positions
    for(set<Alignment*>::iterator iter = set_of_alignments.begin();
        iter != set_of_alignments.end();
        iter++)
    {
        sorted_by_X.insert(make_pair((*iter)->start_x,*iter));
        sorted_by_X.insert(make_pair((*iter)->end_x,  *iter));
    }
    // Assign an optinal alignment to each alignment by updating the list sorted by Y
    for(multimap<int, Alignment*>::iterator iter = sorted_by_X.begin();
        iter != sorted_by_X.end();
        iter++)
    {
        Alignment* tmpX_alignment = iter->second;
        multimap<int, Alignment*>::iterator tmpY, prevY;
        
        if(tmpX_alignment->isStart(iter->first)){
            if(!sorted_by_Y.empty()){
                // Search for a position such that tmpY->end_y <= tmpX_alignment->end_y
                for(tmpY = sorted_by_Y.begin(), prevY = tmpY;
                    tmpY != sorted_by_Y.end();
                    prevY = tmpY, tmpY++)
                {
                    if(prevY->second->end_y <= tmpX_alignment->start_y + MAX_LEN_overlapping_alignments &&
                       tmpY->second->end_y   > tmpX_alignment->start_y + MAX_LEN_overlapping_alignments )
                    {
                        tmpX_alignment->set_predecessor(prevY->second);
                        break;
                    }
                }
                if(prevY->second->end_y <= tmpX_alignment->start_y + MAX_LEN_overlapping_alignments &&
                   tmpY == sorted_by_Y.end())
                {
                    tmpX_alignment->set_predecessor(prevY->second);
                }
            }
        }else{
            
            if(sorted_by_Y.empty()){
                sorted_by_Y.insert(make_pair(tmpX_alignment->end_y, tmpX_alignment));
#ifdef DEBUG_chaining
                cout << "insert\t";
                tmpX_alignment->print_one_TR(print_alignment);
#endif
            }else{
                bool flag = true;   // flag for indicating that tmpX_alignment should be inserted
                for(tmpY = sorted_by_Y.begin(); tmpY != sorted_by_Y.end(); tmpY++)
                {
                    if(tmpY->second->end_y <= tmpX_alignment->end_y &&
                       tmpY->second->score > tmpX_alignment->score)
                    {   // A better alignment was found.
                        flag = false;
                    }
                    if(tmpY->second->end_y > tmpX_alignment->end_y){
                        break;
                    }
                }
                if(flag){
                    sorted_by_Y.insert(make_pair(tmpX_alignment->end_y, tmpX_alignment));
#ifdef DEBUG_chaining
                    cout << "insert\t";
                    tmpX_alignment->print_one_TR(print_alignment);
#endif
                    // Delete unnecessary alignments
                    for(tmpY = sorted_by_Y.begin();
                        tmpY != sorted_by_Y.end();
                        tmpY++)
                    {
                        if(tmpY->second->end_y >= tmpX_alignment->end_y &&
                           tmpY->second->score < tmpX_alignment->score)
                        {
#ifdef DEBUG_chaining
                            cout << "delete\t"; tmpY->second->print_one_TR(print_alignment);
#endif
                            sorted_by_Y.erase(tmpY);
                        }
                    }
                }
            }
        }
    }
    
    // Print the maximum chain
    // Concatenate fragmented alignmnets of an identical tandem repeat
    (sorted_by_Y.rbegin())->second->concatenate_similar_alignments();
    
#ifdef DEBUG_chaining
    (sorted_by_Y.rbegin())->second->print_chain();
#endif
    (sorted_by_Y.rbegin())->second->print_all_TRs(print_alignment);

    // delete all
    for(set<Alignment*>::iterator iter = set_of_alignments.begin();
        iter != set_of_alignments.end(); iter++){
        set_of_alignments.erase(iter);
    }
    for(multimap<int, Alignment*>::iterator iter = sorted_by_X.begin();
        iter != sorted_by_X.end(); iter++){
        sorted_by_X.erase(iter);
    }
    for(multimap<int, Alignment*>::iterator iter = sorted_by_Y.begin();
        iter != sorted_by_Y.end(); iter++){
        sorted_by_Y.erase(iter);
    }
}

