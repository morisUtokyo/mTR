#include "chaining.h"
#include "mTR.h"
#include <stdio.h>
#include <iostream>
#include <set>
#include <map>
#include <string.h>
using namespace std;
//#define DEBUG_chaining

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
    
    void print_one_TR()const{
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
        fflush(stdout);
        
    }
    void print_all_TRs()const{
        if(predecessor != NULL){
            predecessor->print_all_TRs();
        }
        print_one_TR();
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

void search_max(){
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
    max_alignment->print_one_TR();
    
    for(set<Alignment*>::iterator iter = set_of_alignments.begin();
        iter != set_of_alignments.end(); iter++){
        set_of_alignments.erase(iter);
    }
}


void chaining(){
    if(set_of_alignments.empty()){
        return;
    }
    multimap<int, Alignment*> sorted_by_X;
    map<int, Alignment*> sorted_by_Y;
    
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
        map<int, Alignment*>::iterator tmpY, prevY;
        
        if(tmpX_alignment->isStart(iter->first)){
#ifdef DEBUG_chaining
            // Pair of an alignment and its start position
            cout << "\nS " << iter->first << "\t" << tmpX_alignment->start_y << "\t";
            tmpX_alignment->print();
            
            cout << "tmpY";
            if(sorted_by_Y.empty()){
                cout << " is empty\n";
            }else{
                cout << " has\n";
                for(tmpY = sorted_by_Y.begin(); tmpY != sorted_by_Y.end(); tmpY++){
                    cout << tmpY->first << "\t";
                    tmpY->second->print();
                }
            }
#endif
            if(!sorted_by_Y.empty()){
                // Search for a position such that tmpY->end_y <= tmpX_alignment->end_y
                for(tmpY = sorted_by_Y.begin(), prevY = tmpY;
                    tmpY != sorted_by_Y.end();
                    prevY = tmpY, tmpY++){
                    if(prevY->second->end_y <= tmpX_alignment->start_y &&
                       tmpY->second->end_y > tmpX_alignment->start_y)
                    {
                        tmpX_alignment->set_predecessor(prevY->second);
#ifdef DEBUG_chaining
                        cout << "Set predececcor\t"; tmpX_alignment->print();
#endif
                        break;
                    }
                }
                if(prevY->second->end_y <= tmpX_alignment->start_y &&
                   tmpY == sorted_by_Y.end())
                {
                    tmpX_alignment->set_predecessor(prevY->second);
#ifdef DEBUG_chaining
                    cout << "Set predececcor\t"; tmpX_alignment->print();
#endif
                }
            }
        }else{
#ifdef DEBUG_chaining
            // Pair of an alignment and its end position
            cout << "\nE " << iter->first << "\t";
            tmpX_alignment->print();
#endif
            
            if(sorted_by_Y.empty()){
                sorted_by_Y.insert(make_pair(tmpX_alignment->end_y, tmpX_alignment));
#ifdef DEBUG_chaining
                cout << "insert\t" << tmpX_alignment->end_y << "\t";
                tmpX_alignment->print();
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
                }
                if(flag){
                    sorted_by_Y.insert(make_pair(tmpX_alignment->end_y, tmpX_alignment));
#ifdef DEBUG_chaining
                    cout << "insert\t" << tmpX_alignment->end_y << "\t";
                    tmpX_alignment->print();
#endif
                    // Delete unnecessary alignments
                    for(tmpY = sorted_by_Y.begin();
                        tmpY != sorted_by_Y.end();
                        tmpY++){
                        if(tmpY->second->end_y >= tmpX_alignment->end_y &&
                           tmpY->second->score < tmpX_alignment->score)
                        {
#ifdef DEBUG_chaining
                            cout << "delete\t"; tmpY->second->print();
#endif
                            sorted_by_Y.erase(tmpY);
                        }
                    }
                }
            }
#ifdef DEBUG_chaining
            cout << "tmpY";
            if(sorted_by_Y.empty()){
                cout << " is empty\n";
            }else{
                cout << " has\n";
                for(tmpY = sorted_by_Y.begin(); tmpY != sorted_by_Y.end(); tmpY++){
                    cout << tmpY->first << "\t";
                    tmpY->second->print();
                }
            }
#endif
        }
    }
    
    // Print the maximum chain
#ifdef DEBUG_chaining
    (sorted_by_Y.rbegin())->second->print_chain();
#endif
    (sorted_by_Y.rbegin())->second->print_all_TRs();
    
    // delete all
    for(set<Alignment*>::iterator iter = set_of_alignments.begin();
        iter != set_of_alignments.end(); iter++){
        set_of_alignments.erase(iter);
    }
    for(multimap<int, Alignment*>::iterator iter = sorted_by_X.begin();
        iter != sorted_by_X.end(); iter++){
        sorted_by_X.erase(iter);
    }
    for(map<int, Alignment*>::iterator iter = sorted_by_Y.begin();
        iter != sorted_by_Y.end(); iter++){
        sorted_by_Y.erase(iter);
    }
}

