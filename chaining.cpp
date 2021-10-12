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

    char*   readID;
    //char    readID[BLK];
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
    int     match_gain;
    int     mismatch_penalty;
    int     indel_penalty;
    char    *string;
    int     *string_score;

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
              int   a_match_gain,
              int   a_mismatch_penalty,
              int   a_indel_penalty,
              char* a_string,
              int*  a_string_score)
    {
        readID      = new char[BLK];
        string      = new char[BLK];
        string_score= new int[MAX_PERIOD];
        
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
        Kmer            = a_Kmer;
        match_gain      = a_match_gain;
        mismatch_penalty= a_mismatch_penalty;
        indel_penalty   = a_indel_penalty;
        strcpy( string, a_string );
        for(int i=0; i<rep_period; i++){
            string_score[i] = a_string_score[i]; }
        
	}
    
    ~Alignment(){
        delete [] readID;
        delete [] string;
        delete [] string_score;
    }
    
    void print_one_TR(int print_alignment)const{
        
        printf(
           "%.50s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%s\n",
           readID,
           inputLen,
           rep_start+1, // 1-origin
           rep_end+1,   // 1-origin
           repeat_len,
           rep_period,
           Num_freq_unit,
           Num_matches,
           (float)Num_matches/repeat_len,
           Num_mismatches,
           Num_insertions,
           Num_deletions,
           string
        );
        
#ifdef DEBUG_unit_score
    printf("\n\t%s", string);
    for(int k=Kmer; Kmer<=k; k--){
        printf("\n%i\t", k);
        print_freq(rep_start, rep_end, rep_period, const_cast<char*>(string), inputLen, k);
    }
    printf("\n");
#endif


        
#ifdef Print_overlapping_event
        if(predecessor != NULL){
            if(predecessor->rep_end  > rep_start){
                printf("------------ overlapping ----------------\n");
            }
        }
#endif
        
        if(print_alignment == 1){
            printf("\n");
            pretty_print_alignment( const_cast<char*>(string), rep_period, rep_start, rep_end, match_gain, mismatch_penalty, indel_penalty);
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
                      int   match_gain,
                      int   mismatch_penalty,
                      int   indel_penalty,
                      char* string,
                      int*  string_score)
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
                            match_gain,
                            mismatch_penalty,
                            indel_penalty,
                            string,
                            string_score));
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
        if((*iter)->start_x + MAX_LEN_overlapping <=  (*iter)->end_x ){
            sorted_by_X.insert(make_pair( (*iter)->start_x,*iter) );
            sorted_by_X.insert(make_pair( ((*iter)->end_x - MAX_LEN_overlapping),  *iter) );
        }
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
                    if(prevY->second->end_y <= tmpX_alignment->start_y + MAX_LEN_overlapping &&
                       tmpY->second->end_y   > tmpX_alignment->start_y + MAX_LEN_overlapping )
                    {
                        tmpX_alignment->set_predecessor(prevY->second);
                        break;
                    }
                }
                if(prevY->second->end_y <= tmpX_alignment->start_y + MAX_LEN_overlapping &&
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
    //(sorted_by_Y.rbegin())->second->concatenate_similar_alignments();
    
#ifdef DEBUG_chaining
    (sorted_by_Y.rbegin())->second->print_chain();
#endif
    (sorted_by_Y.rbegin())->second->print_all_TRs(print_alignment);

    
    // delete all
    /*
    set_of_alignments.clear();
    sorted_by_X.clear();
    sorted_by_Y.clear();
    */
    for(multimap<int, Alignment*>::iterator iter = sorted_by_X.begin();
        iter != sorted_by_X.end(); iter++){
        sorted_by_X.erase(iter++);
    }
    for(multimap<int, Alignment*>::iterator iter = sorted_by_Y.begin();
        iter != sorted_by_Y.end(); iter++){
        sorted_by_Y.erase(iter++);
    }
    for(set<Alignment*>::iterator iter = set_of_alignments.begin();
        iter != set_of_alignments.end(); iter++){
        // delete all elements of Alignment
        Alignment *tmp = (*iter);
        delete tmp;
        set_of_alignments.erase(iter++);
    }
}

