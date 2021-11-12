
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <iterator>

// All headers included with quote marks are local files. All headers with <> are part of C++ standard library.
#include "sealcensusFullCpp.h"

using namespace std;



/*   
This has the main function for executing the sealcensus estimator. All files needed for compiling are in the same folder. Besides the 7 header files (extension .h), there is utilities.cpp, and one file with sample data (SampleData.csv). There is a folder within named files where the entire chain of parameter estimates for every season are stored during program execution.  

Compiling requires the standard linux C++ compiler, gcc, updated recently (some early versions will fail). This linux command will compile the executable: 
g++ -Wall sealcensusFull.cpp -o sealcensus.exe

Then to execute, 5 arguments must be appended to the sealcensus.exe: the name of the data file, the input parameter file, the number of steps to execute, the number of steps to discard as burn-in, and an integer to indicate how often progress should be printed to the screen. For example:
./sealcensus.exe SampleData.csv inputfile.txt 6000 1000 999
*/

/* The main function accepts command-line parameters, then invokes a SealCensus class, passing the needed arguments. The class is defined in sealcensusFullCpp.h */
int main(int argc, char *argv[])
{
  string infile=argv[1], startfile=argv[2];
  int burn=atoi(argv[4]), steps=atoi(argv[3]), show=atoi(argv[5]);
  IntVector stepping{steps,burn,show};  
  
  time_t rightnow = time(NULL);
  printf("Start reading files %s and %s at %s", infile.c_str(), startfile.c_str(), ctime(&rightnow));
  
  SealCensus C(infile,startfile,stepping);
  C.FullRun();
  
  rightnow = time(NULL);
  printf("Finished %d steps at %s", steps, ctime(&rightnow));
}
