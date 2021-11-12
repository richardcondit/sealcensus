#ifndef READ_METROP_H
#define READ_METROP_H 1

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <iterator>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <random>
#include "/home/condit/programs/utilities/statistics.h"
#include "/home/condit/programs/utilities/likelihood.h"
#include "/home/condit/programs/utilities/util.h"
#include "/home/condit/programs/utilities/utilities.cpp"
#include "/home/condit/programs/utilities/Array2D.h"

using namespace std;

/*   
C++ class whose sole purpose is reading an ascii table of successive Metropolis parameter estimates. The input includes a file name holding columns of parameters, the number of steps to ignore as burn-in. The number of columns and rows is arbitrary, worked out by reading the file. Results are written to a new file whose name is derived from the input file, which is assumed to have extension csv.
 
cd ~/elephantseal/census.analysis/cpp
g++ -Wall readMetropParam.cpp -o readParam.exe 2>&1 | more
g++ -Wall readMetropParam.cpp -o readParam.exe
./readParam.exe Year2018 500
*/


class ReadParam
   {	
	protected:
        int burnin, postburn=0, Ncol, N=0;
        Array2D<NumVector> fullparam, sortparam;
        unordered_map<string,NumVector> summary;

	public:
        ReadParam(string, int); 
        void PrintResult(void);
        void SortParam(void);
        unordered_map<string,NumVector> CalcSummary(void);
   };

// Initialize the class, requiring the filename with parameters and the number of burn-in steps. Post-burn-in rows are read into a 2D array of double (ArrayD<NumVector>).
ReadParam::ReadParam(string filename, int B)
   {
	ifstream infile(filename);	
	burnin=B;
	string oneline;
	getline(infile,oneline);   // Header row
	StrVector HCells;
	TokenizeStr(oneline, HCells, "\t");
	Ncol=(int)(HCells.size());
	fullparam.SetSize(0,Ncol);
	  
	while(getline(infile,oneline))
	 {
		StrVector Cells;
		TokenizeStr(oneline, Cells, "\t");
		NumVector oneset;
		for(int j=0;j<Ncol;j++) oneset.push_back(stod(Cells[j]));
		if(N>burnin) fullparam.AddOneRow(oneset);
		N++;
	 }

	postburn=fullparam.rows();
	if(postburn>0) SortParam();	
   }

   
// Because percentiles require sorted, every column of fullparam is sorted into a new Array2D.
void ReadParam::SortParam()
{
	sortparam.SetSize(postburn,Ncol);
    for(int j=0;j<Ncol;j++) 
     {
		 NumVector onecol=fullparam.GetOneCol(j);
		 sort(onecol.begin(),onecol.end());
		 sortparam.SetOneCol(j,onecol);
	 }
   	
}

// Calculate mean, median, and 95th percentiles of every column of fullparam. 
unordered_map<string,NumVector> ReadParam::CalcSummary()
{
	if(postburn==0) return(summary);
	
	NumVector lower, upper, mean, median;
	int percent5=(int)(round(postburn)/20);
	int middle=(int)(round(postburn)/2);
	
    for(int j=0;j<Ncol;j++)	
     {
		 NumVector onecol=sortparam.GetOneCol(j);
		 lower.push_back(*(onecol.begin()+percent5));
		 upper.push_back(*(onecol.end()-percent5));
		 median.push_back(*(onecol.begin()+middle));
		 mean.push_back(Mean(onecol));
	 }
	 
	summary.insert(make_pair("Mean",mean));
	summary.insert(make_pair("Median",median));
	summary.insert(make_pair("Lower",lower));
	summary.insert(make_pair("Upper",upper));	
    return(summary);	
}



/* For testing purposes, a standalone version. Must be commented to include this file elsewhere. 
int main(int argc, char *argv[])
{
  string infile=argv[1];
  int burn=atoi(argv[2]);

  ReadParam param(infile,burn);
}
*/
 

#endif
