#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <strstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

class reads{
	
public:
	int ID;
	int length, start, end, rep_length, unit_length, unit_num;
	string seq;
	int match, mismatch, ins, del, kmer, method;
	double match_ratio;
};

vector<string> read_repeat(const char* filename)
{
	string line;

	ifstream fin; 

	fin.open(filename); 
	if (!fin){ cout << "error : file not found" << endl; exit(1); } // if fin==0  file not exist

	vector<string> rep_row;
	
	while (getline(fin, line)){
		rep_row.push_back(line);
	}

	return rep_row;
}

void read_out(const char* filename, vector<string> rep_row)
{
	string line;

	ifstream fin; 

	fin.open(filename); 
	if (!fin){ cout << "error : file not found" << endl; exit(1); } // if fin==0  file not exist
	
	int count = 0;
	
	vector<int> perfect_row;
	
	for(int i=0; i<(int)rep_row.size(); i++)
	{
		perfect_row.push_back(0);
	}
	
	while (getline(fin, line)){

		count++;
		
		if (line.empty())
		{
			break;
		}
		
		replace(line.begin(), line.end(), '	', ' '); // tab to space
		replace(line.begin(), line.end(), ',', ' '); // tab to space
		replace(line.begin(), line.end(), ')', ' '); // tab to space
		
		reads n_reads;
		
		stringstream buf(line);
		
		buf >> n_reads.ID >> n_reads.length >> n_reads.start >> n_reads.end >> n_reads.rep_length >> n_reads.unit_length >> n_reads.unit_num >> n_reads.match >> n_reads.match_ratio >> n_reads.mismatch >> n_reads.ins >> n_reads.del >>  n_reads.kmer >> n_reads.method >> n_reads.seq;

		string double_sequence;
		double_sequence = rep_row[n_reads.ID] + rep_row[n_reads.ID];
		
		int max = 0;
		
		for(int i=0; i<(int)double_sequence.size(); i++)
		{
			int match = 0;
			
			for(int j=0; j<(int)n_reads.seq.size(); j++)
			{
				if(i+j < (int)double_sequence.size())
				{
					if(double_sequence[i + j] == n_reads.seq[j])
					match++;
				}
			}
			if(max < match)
			{
				max = match;
			}
		}
		
		if(max == (int)rep_row[n_reads.ID].size())
		{
			perfect_row[n_reads.ID] = 1;
		}
	}
	
	int count_p_all = 0;
	
	for(int i=0; i<(int)perfect_row.size(); i++)
	{
		if(perfect_row[i] == 1)
		{
			count_p_all++;
		}
	}
	
	cout << count_p_all << endl;
}

int main(int argc, char** argv)
{
	vector<string> rep_row = read_repeat(argv[2]);
	read_out(argv[1], rep_row);
}

