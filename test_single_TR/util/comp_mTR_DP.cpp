//
//  count_match.cpp
//
//
//  Created by Kazuki Ichikawa
//
#include <iostream>
#include <fstream>
#include <string>
#include <strstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>

using namespace std;

double threshold = 0.01;

class mTR{	
public:
	int ID;
	int length, start, end, rep_length, unit_length, unit_num;
	string seq;
	int match, mismatch, ins, del;
	double match_ratio;
};

string rev_atcg(char s)
{
	if(s == 'A')
	{
		return "T";
	}
	else if(s == 'T')
	{
		return "A";
	}
	else if(s == 'C')
	{
		return "G";
	}
	else if(s == 'G')
	{
		return "C";
	}
	else
	{
		return "N";
	}
}

string rev_str(string s)
{
	string tmp = "";
	for(int i=0; i<(int)s.size(); i++)
	{
		tmp += rev_atcg(s[(int)s.size() - 1 - i]);
	}
	return tmp;
}

double cal_DP(string a, string b)
{
	int match = 1;
	int miss = -1;
	int gap = -1;

	vector< vector <int> > matrix( ((int)b.length() + 1), vector<int>( ((int)a.length() + 1) , -999999) );

	matrix[0][0] = 0;

	for(int i=0; i< (int)matrix.size(); i++)
	{
		matrix[i][0] = 0;
	}

	for(int j=1; j < (int)matrix[0].size(); j++)
	{
		for(int i=0; i < (int)matrix.size(); i++)
		{
			if(i != 0)
			{
				if(matrix[i][j] < matrix[i - 1][j] + gap)
				{
					matrix[i][j] = matrix[i - 1][j] + gap;
				}

				if(j != 0)
				{
					if(a[j - 1] == b[i - 1])
					{
						if(matrix[i][j] < matrix[i - 1][j - 1] + match)
						{
							matrix[i][j] = matrix[i - 1][j - 1] + match;
						}
					}
					else
					{
						if(matrix[i][j] < matrix[i - 1][j - 1] + miss)
						{
							matrix[i][j] = matrix[i - 1][j - 1] + miss;
						}						
					}
				}
			}
			else
			{
				int last = (int)b.length() - 1;

				if(j != 0)
				{
					if(a[j - 1] == b[last])
					{
						if(matrix[i][j] < matrix[(int)matrix.size() - 2][j - 1] + match)
						{
							matrix[i][j] = matrix[(int)matrix.size() - 2][j - 1] + match;
						}
					}
					else
					{
						if(matrix[i][j] < matrix[(int)matrix.size() - 2][j - 1] + miss)
						{
							matrix[i][j] = matrix[(int)matrix.size() - 2][j - 1] + miss;
						}						
					}
				}
			}

			if(j != 0)
			{
				if(matrix[i][j] < matrix[i][j - 1] + gap)
				{
					matrix[i][j] = matrix[i][j - 1] + gap;
				}
			}
		}
	}

	//print_matrix(matrix);

	string r_a = "";
	string r_b = "";

	int x = (int)matrix.size() - 1;
	int y = (int)matrix[0].size() - 1;

	int max = matrix[x][y];

	for(int i=0; i < (int)matrix.size(); i++)
	{
		if(matrix[i][y] > max)
		{
			max = matrix[i][y];
			x = i;
		}
	}

	int start_x = x;
	int match_num = 0;

	while(1)
	{
		bool update = false;

		if(x == 0 && y > 0)
		{
			int last = (int)b.length() - 1;

			if((a[y - 1] == b[last]) && (matrix[x][y] - match == matrix[(int)matrix.size() - 2][y - 1]))
			{
				r_a += a[y - 1];
				r_b += b[last];

				x = (int)matrix.size() - 2;
				y--;
				update = true;
				match_num++;
			}
			else if((a[y - 1] != b[last]) && (matrix[x][y] - miss == matrix[(int)matrix.size() - 2][y - 1]))
			{
				r_a += a[y - 1];
				r_b += b[last];

				x = (int)matrix.size() - 2;
				y--;
				update = true;				
			}
		}

		if(x > 0 && y > 0 && !update)
		{
			if((a[y - 1] == b[x - 1]) && (matrix[x][y] - match == matrix[x - 1][y - 1]))
			{
				r_a += a[y - 1];
				r_b += b[x - 1];

				x--; y--;
				update = true;
				match_num++;
			}
			else if( (a[y - 1] != b[x - 1]) &&  (matrix[x][y] - miss == matrix[x - 1][y - 1]))
			{
				r_a += a[y - 1];
				r_b += b[x - 1];

				x--; y--;
				update = true;
			}
		}

		if(x > 0 && !update)
		{
			if(matrix[x][y] - gap == matrix[x - 1][y])
			{
				r_a += "-";
				r_b += b[x - 1];

				x--;
				update = true;
			}
		}

		if(y > 0 && !update)
		{
			if(matrix[x][y] - gap == matrix[x][y - 1])
			{
				r_a += a[y - 1];
				r_b += "-";

				y--;
			}
		}

		if(y == 0)
		{
			break;
		}
	}

	reverse(r_a.begin(), r_a.end());
	reverse(r_b.begin(), r_b.end());

	/*cout << r_a << endl;

	for(int i=0; i<(int)r_a.length(); i++)
	{
		if(r_a[i] == r_b[i])
		{
			cout << "|";
		}
		else
		{
			cout << " ";
		}
	}
	cout << endl;

	cout << r_b << endl;

	cout << match_num << endl;

	int check;
	cin >> check;*/

	double match_ratio = (double)match_num/(double)r_a.length();
	return match_ratio;
}

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

void read_mTR(const char* filename, const char* out_name, vector<string> rep_row)
{
	string line;

	ifstream fin; 

	fin.open(filename); 
	if (!fin){ cout << "error : file not found" << endl; exit(1); } // if fin==0  file not exist

	ofstream ofs;
	ofs.open(out_name);

	int n = 0;

	int dif_n = (int)(rep_row[0].length() * threshold);

	while (getline(fin, line)){
		
		if (line.empty())
		{
			break;
		}

		double match_ratio = 0;

		mTR n_mTR;
		stringstream buf(line);

		buf >> n_mTR.ID >> n_mTR.length >> n_mTR.start >> n_mTR.end >> n_mTR.rep_length >> n_mTR.unit_length >> n_mTR.unit_num >> n_mTR.match >> n_mTR.match_ratio >> n_mTR.mismatch >> n_mTR.ins >> n_mTR.del >> n_mTR.seq;

		string a, b;

		if(rep_row[n_mTR.ID].length() >= n_mTR.seq.length())
		{
			a = rep_row[n_mTR.ID]; b = n_mTR.seq;
		}
		else
		{
			a = n_mTR.seq; b = rep_row[n_mTR.ID];
		}

		match_ratio = cal_DP(a, b);

		/*if((double)(b.length() *(1+threshold)) >= (double)a.length())
		{
			cout << n << endl;
			match_ratio = cal_DP(a, b);

			string rev_s = rev_str(a);
			double tmp = cal_DP(rev_s, b);

			if(tmp > match_ratio)
			{
				match_ratio = tmp;
			}
		}*/

		ofs << match_ratio << endl;
		//ofs << match_ratio << "	" << line << endl;

		/*if(match_ratio >= 0.8)
		{
			ofs << line << endl;
		}*/

	}
}

int main(int argc, char** argv)
{
	vector<string> rep_row = read_repeat(argv[1]);
	read_mTR(argv[2], argv[3], rep_row);
}