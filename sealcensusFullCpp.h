#ifndef SEALCENSUS_FULL_RCPP_H
#define SEALCENSUS_FULL_RCPP_H 1

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <iterator>
#include <map>
#include <set>
#include <algorithm>
#include <random>

// All headers included with quote marks are local files. All headers with <> are part of C++ standard library.
#include "randomgenerator.h"
#include "statistics.h"
#include "utilities.cpp"
#include "Array2D.h"
#include "readMetropParam.h"

using namespace std;

/*   
This file defines the SealCensus class, which executes the program to estimate population size from multiple censuses on days throughout one or more seasons. 
*/

typedef map<int,Array2D<IntVector>> Int2IntArray;  
typedef map<int,NumVector> Int2NumVect;
typedef map<int,string> Int2Str;


class SealCensus
   {	
	protected:
		RandomGenerator rand;
	    int counter=0, sample=0, Npar; // Step counter, number of observations, number of parameters
 	    string infilename, paramfile;  // infilename names file with data, paramfile has priors and day range
 	    string hyperfile, hyperfileSD; // One file for chain of hyper-parameters, and one for hyperSD parameters
 	    Int2Str yearfile;              // One file name for each year, each file with full chain parameters; indexed to integer year
	    IntVector stepping, dayrange;  // dayrange read from paramfile -- dayrange = start and end day; stepping is input with initialization: steps, burnin, show
	    double target=0.25, adjexp, adjust=1.01, fullLikelihood, hyperLikelihood;
	    Array2D<IntVector> obs;        // Input observations, 3 columns of integers, year, day, count

		set<int> allyear;              // Unique years stored as a set, then used as index for Xlst
		Int2IntArray Xlst;             // Map format for data
		NumVector startparam{0,0,0,0,0,0}, scale, hyper, hyperSD, scaleHyper, scaleHyperSD;
		Int2NumVect fullparam;         // One set of parameters per year, written to file so not stored

		// Hard-coded priors, but can be replaced inputfile.txt. If that file has a zero in first entry of second row, these default are used.
		double tenurePriorMean=31.06, tenurePriorSD=0.27, cvTenurePriorMean=0.1212, cvTenurePriorSD=0.0062;	    
	    ofstream hyperSDstream, hyperstream, yearstream, output; // File streams for parameter chains. The yearstream handles the entire vector of yearfile.
    	bool showit=0;

	public:
	    /* The class constructor accepts 3 arguments and assigns them to class variables:
	    filename is data file, S is size 3, steps to run with burnin and show, range is date range, set prior is 4 integers for priors */
        SealCensus(string datafile, string inputfile, IntVector S): infilename(datafile), paramfile(inputfile), stepping(S) {}; 
		/* The class function FullRun controls execution, calling all of the other class functions */
		void FullRun(void);

		void Initialize(void);
        void InitializePriors(void);
		void FillMap(void);
		void SetStartParam(void);
		void updateOneParam(int,int);
		void updateOneHyper(int);
		void updateOneHyperSD(int);
		void UpdateParam(void);
		double llikeSealCensus(Array2D<IntVector>,NumVector);
		double hyperllike(double, int);
        Array2D<NumVector> PredictSealCensus(NumVector);
        NumVector censusModel(IntVector,NumVector);
        void StartPrint(void);
        void PrintHyper(void);
        void PrintYears(void);
        void ReturnSummary(void);
        void PrintOneResult(unordered_map<string,NumVector>,string,string);
   };

/* The class function FullRun controls execution by calling all the other class functions. After class construction, it is the only class function called by main. 
 * The first 5 subroutines set up the data, parameters, and open files for printing the chain of parameter estimates.
 * Then a loop is invoked with the requested number of steps for the chain (the class variable stepping[0]). Within the loop, parameters are update (UpdateParam), then printed to the growing files by PrintHyper and PrintYears.
 * The cout line outputs the current value of hyper-parameters to the screen ever stepping[2] steps.
 * After finishing the loop, the two files with hyper-parameters are closed, and the final subroutine called to calculate all summary statistics (ReturnSummary).
 */
void SealCensus::FullRun()
{
	Initialize();
	InitializePriors();
	FillMap();
	SetStartParam();
    StartPrint();
    
	for(int i=0;i<stepping[0];i++) 
	 {
	   UpdateParam();
	   PrintHyper();
	   PrintYears();
	   
       counter++;
	   if(counter%stepping[2]==0) showit=1;  // Debug printing to screen
       else showit=0;

       if(showit) { cout << "Step " << counter << ", hyper: "; Test(hyper); }
     } 
   
   hyperstream.close();    // These two streams remain open throughout the run, must be closed at the end.
   hyperSDstream.close();
   ReturnSummary();
}

	

// The class function Initialize reads the file with census data. The name of the file is submitted at execution time on the command line. The data are stored in a class variable obs, a 2D array of integers. The number of rows read is stored in the class variable sample. There are 3 columns.
void SealCensus::Initialize()
{
	ifstream infile(infilename), startfile(paramfile);
	string oneline;
	getline(infile,oneline);   // Header row
	IntVector year, day, count;
    
	while(getline(infile,oneline))
	 {
		StrVector Cells;
		TokenizeStr(oneline, Cells, "\t");
		int y=stoi(Cells[0]);
		int d=stoi(Cells[1]);
		int c=stoi(Cells[2]);
		
		year.push_back(y);
		day.push_back(d);
		count.push_back(c);
		sample++;
	 }
	
   obs.SetSize(sample,3);   // Array of sample rows, 3 cols
   obs.SetOneCol(0,year);
   obs.SetOneCol(1,day);
   obs.SetOneCol(2,count);
}

// The class function InitializePriors reads a file having 6 input parameters on 2 rows. The first row has the date range, and the second has the 4 necessary priors: tenurePriorMean, tenurePriorSD, cvTenurePriorMean, cvTenurePriorSD;	    
void SealCensus::InitializePriors()
{
	ifstream startfile(paramfile);
	string oneline;
	StrVector Row1, Row2;

	getline(startfile,oneline);   
	TokenizeStr(oneline, Row1, "\t");	
	dayrange.push_back(stoi(Row1[0]));
	dayrange.push_back(stoi(Row1[1]));
	
	getline(startfile,oneline); 
	TokenizeStr(oneline, Row2, "\t");	
	double firstentry=stod(Row2[0]);
	if(firstentry>0)
	 {
      tenurePriorMean=stod(Row2[0]);
      tenurePriorSD=stod(Row2[1]);
      cvTenurePriorMean=stod(Row2[2]);
      cvTenurePriorSD=stod(Row2[3]);
     }
}

// The class function SetStartParam calculates the initial parameter values for the Metropolis chains. There must be a set of parameters for every season, plus a set of hyper-parameters defining variation across seasons. This is run after FillMap fills the variable Xlst with all the observations. The starting parameters are calculated using the mean of the maximum count per year, its variance, and the mean and variance of the day each year. These are not good estimates, but they are in the right ballpark.
void SealCensus::SetStartParam()
{
	double SDscale;
	NumVector maxpop, meanday, sdday;
	for(auto it=allyear.begin();it!=allyear.end();it++) 
	 {
	  meanday.push_back(Mean(Xlst[*it].GetOneCol(1)));
	  sdday.push_back(sqrt(Var(Xlst[*it].GetOneCol(1))));
	  maxpop.push_back(Max(Xlst[*it].GetOneCol(2)));
	 }
	startparam[0]=Mean(maxpop);
	startparam[1]=Mean(meanday);
	startparam[2]=Mean(sdday);	
	startparam[3]=0;
	
	startparam[4]=tenurePriorMean;
	startparam[5]=cvTenurePriorMean;
	Npar=(int)startparam.size();

    for(auto it=allyear.begin();it!=allyear.end();it++) fullparam.insert(make_pair(*it,startparam));
	
	for(auto it=startparam.begin();it!=startparam.end();it++) 
	 { 
	   if(*it>0) SDscale=(*it)/10;
	   else SDscale=0.1;
	   scale.push_back(SDscale); 

	   hyper.push_back(*it);
	   hyperSD.push_back(SDscale); 
	   scaleHyper.push_back(SDscale);
       scaleHyperSD.push_back(SDscale);
	 }
	 
   adjexp=(1-target)/target;
}

// The input data in the variable obs are converted to a map object stored in the class variable Xlst. The index of Xlst is the seasons (ie, years). This allows data for a single year to be accessed easily. The class function FillMap fills Xlst from obs. It is only needed once, at the outset.  
void SealCensus::FillMap()
{
   Array2D<IntVector> empty(0,obs.cols());

   for(int i=0;i<sample;i++)
    {
	  IntVector onerow=obs[i];
	  allyear.insert(onerow[0]);
	}

   for(auto it=allyear.begin();it!=allyear.end();it++) Xlst.insert(make_pair(*it,empty));
    	  	
   for(int i=0;i<sample;i++)
    {
	  int yr=obs[i][0];
	  Xlst[yr].AddOneRow(obs[i]);
	}
	  	
}

// The class function UpdateParam is a cover, simply calling the three other update functions. It is executed once each step of the Metropolis chain. 
void SealCensus::UpdateParam()
{
   for(auto it=allyear.begin();it!=allyear.end();it++) for(int j=0;j<Npar;j++) updateOneParam(*it,j);
   for(int j=1;j<Npar;j++) updateOneHyper(j);  
   for(int j=1;j<Npar;j++) updateOneHyperSD(j);  
}


/* The 3 functions to update parameters are the key of the Metropolis process. Each is called once for every step of the execution. 
 * First, the likelihood functions are called to calculate the log-likelihood of the current set of parameters, one year at a time, but including hyper-parameters. 
 * Then just one parameter is updated to a new value drawn at random (the function rand.norm from randomgenerator.h). In updateOneParam, it is a single parameter for one of the years. Each year is treated separately from all the others, but they all depend on hyper-parameters.
 * A new log-likelihood is calculated using all the parameters, using the current (unchanged) value for all except the one that was updated.
 * The Metropolis algorithm is the decision about whether to accept the new parameter, depending on the difference between the original and new log-likelihoods.
 * Successive parameters of the Metropolis chain are not stored, they are simply printed to expanding files. Thus, if accepted, the one parameter is simply substituted with the updated value.
 */
void SealCensus::updateOneParam(int y, int w)
{
	double newval, origlike, censlike, newlike, likeratio;
	IntVector days=Xlst[y].GetOneCol(1);
    NumVector par=fullparam[y];
    double h=hyperllike(par[w],w);
    if(w==0) h=0;
    censlike=llikeSealCensus(Xlst[y],par);
    origlike=censlike+h;
    // if(showit && y==1982) printf("Was %6.3lf with likelihoods: census=%6.3lf, hyper=%6.3lf, ", par[w], censlike, h);

	newval=rand.norm(1,par[w],scale[w])[0];
	par[w]=newval;
    if(w>0) h=hyperllike(par[w],w);
    else h=0;
    censlike=llikeSealCensus(Xlst[y],par);
    newlike=censlike+h;
	likeratio=exp(newlike-origlike);
	if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    // if(showit && y==1982) printf("try %6.3lf with likelihoods: census=%6.3lf, hyper=%6.3lf\n", par[w], censlike, h);
    
	NumVector r=rand.unif(1,0.0,1.0);
	if(r[0]<likeratio)                                     // Accept
	 {
	    fullparam[y][w]=newval; 
	    scale[w]*=pow(adjust,adjexp);
	    fullLikelihood=newlike;
     }
	else 
	 {
	    scale[w]*=(1/adjust);                              // Reject		
	    fullLikelihood=origlike;
	 }   
}

// The updateOneHyper function is for updating one of the hyper-means. It requires the parameter estimates for every year.
void SealCensus::updateOneHyper(int w)
{
	double newval, origlike, newlike, likeratio;
	NumVector current;
    for(auto it=allyear.begin();it!=allyear.end();it++) current.push_back(fullparam[*it][w]);
    origlike=Sum(pdfNorm(current,hyper[w],hyperSD[w],1));
    newval=rand.norm(1,hyper[w],scaleHyper[w])[0];
    newlike=Sum(pdfNorm(current,newval,hyperSD[w],1));
	likeratio=exp(newlike-origlike);
	if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    
	NumVector r=rand.unif(1,0.0,1.0);
	if(r[0]<likeratio)                               // Accept
	 {
	    hyper[w]=newval; 
	    scaleHyper[w]*=pow(adjust,adjexp);
	 }
	else scaleHyper[w]*=(1/adjust);                  // Reject		
}


// The updateOneHyperSD function is for updating one of the hyper-standard-deviations. It requires the parameter estimates for every year.
void SealCensus::updateOneHyperSD(int w)
{
	double newval, origlike, newlike, likeratio;
	NumVector current;
    for(auto it=allyear.begin();it!=allyear.end();it++) current.push_back(fullparam[*it][w]);
    // if(showit) Test(current);
    origlike=Sum(pdfNorm(current,hyper[w],hyperSD[w],1));
    newval=rand.norm(1,hyperSD[w],scaleHyperSD[w])[0];
    newlike=Sum(pdfNorm(current,hyper[w],newval,1));
	likeratio=exp(newlike-origlike);
	if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    
	NumVector r=rand.unif(1,0.0,1.0);
	if(r[0]<likeratio)                                // Accept
	 {
	    hyperSD[w]=newval; 
	    scaleHyperSD[w]*=pow(adjust,adjexp);
	    hyperLikelihood=newlike;
	 }
	else 
	 {
		scaleHyperSD[w]*=(1/adjust);                  // Reject		
		hyperLikelihood=origlike;
	 }
}

// The class function hyperllike is a cover to call pdfNorm (from likelihood.h) to calculate the log-likelihood of a single hyperparameter, given the hyper-mean and hyper-standard-deviation.
double SealCensus::hyperllike(double newpar, int w)
{
	double llike=pdfNorm(newpar,hyper[w],hyperSD[w],1);
	return llike;
}


// The class function llikeSealCensus calculates the log-likelihood of the observed counts in a single year (in the submitted Array x), given a set of model parameters (in the NumVector par). It calls the class function censusModel to find predicted censuses on each day with observations. Then it uses poissLogLike (in likelihood.h) to find the Poisson probability of the observed counts given predictions. Since there are prior probabilities assigned to two of the parameters (par[4] and par[5]), pdfNorm must be called for those and included in the total likelihood.
double SealCensus::llikeSealCensus(Array2D<IntVector> x, NumVector par)
{
	 NumVector predlambda=censusModel(x.GetOneCol(1),par);
     IntVector count=x.GetOneCol(2);
     
     if(par[3]<(-.5)) return((-1)*INFINITY);
     if(par[3]>(.35)) return((-1)*INFINITY);
	  
     double tenurePriorllike=pdfNorm(par[4],tenurePriorMean,tenurePriorSD,1);
     double cvPriorllike=pdfNorm(par[5],cvTenurePriorMean,cvTenurePriorSD,1);
   
     double obsllike=poissLogLike(&count,&predlambda);
     double total=tenurePriorllike+cvPriorllike+obsllike;
     
     return(total);
}



// There are two class functions for calculating the predicted number of animals on each day given the model parameters. The first, PredictSealCensus, does the calculation for every day in the entire range of days for a season (in the class variable dayrange). The second function subset from that result the predicted number on the days with observations, since those are the only ones needed for comparing to observations.
Array2D<NumVector> SealCensus::PredictSealCensus(NumVector param)
{
 double N=param[0];           // Population size in one season
 double meanArrival=param[1];
 double sdArrival=param[2];
 double corrstay=param[3];
 double day45stay=param[4];
 double cvTenure=param[5];

 double meanDeparture=day45stay+meanArrival;
 double sdTenure=cvTenure*day45stay;
 double varDeparture=pow(1.0+corrstay,2)*pow(sdArrival,2)+pow(sdTenure,2);
 double sdDeparture=sqrt(varDeparture);
 
 double Ntoday, CumulativeArr=0, CumulativeDep=0;
 NumVector day, CumArrived, NArrival, Arrival, CumDeparted, NDeparture, Departure, PredCensus;
 
 for(int i=dayrange[0];i<=dayrange[1];i++) day.push_back((double)i); // Observation days as double for calculations
 Arrival=pdfNorm(day,meanArrival,sdArrival,0);
 Departure=pdfNorm(day,meanDeparture,sdDeparture,0);
 int len=dayrange[1]-dayrange[0]+1;
 
 for(auto i=0;i<len;i++) 
  {
	Ntoday=N*Arrival[i];
	NArrival.push_back(Ntoday);
	CumulativeArr+=Ntoday;
	CumArrived.push_back(CumulativeArr);
	
    Ntoday=N*Departure[i];
    NDeparture.push_back(Ntoday);
    CumulativeDep+=Ntoday;
    CumDeparted.push_back(CumulativeDep);
    
    PredCensus.push_back(CumulativeArr-CumulativeDep);
  }
   
 Array2D<NumVector> final(PredCensus.size(),6);
 final.SetOneCol(0,day);
 final.SetOneCol(1,PredCensus);
 final.SetOneCol(2,NArrival);
 final.SetOneCol(3,NDeparture);
 final.SetOneCol(4,CumArrived);
 final.SetOneCol(5,CumDeparted);
 
 return(final);
}


// The function that subset predicted population on the days with observations. 
NumVector SealCensus::censusModel(IntVector days, NumVector param)
{
	Array2D<NumVector> fullCensus=PredictSealCensus(param);
	NumVector pred, fullpred=fullCensus.GetOneCol(1), doubleDay=fullCensus.GetOneCol(0);
	IntVector predDay(doubleDay.begin(),doubleDay.end());  // Create an IntVector of every day of the season.
	IntVector::iterator itday, itfull;
	
	for(itday=days.begin();itday!=days.end();itday++)
	  {
		 itfull=find(predDay.begin(),predDay.end(),(*itday));  // Find the prediction for each day with an observation.
		 int steps=distance(predDay.begin(),itfull);
		 pred.push_back(fullpred[steps]);
	  }
	
	return(pred);
}


// The class function StartPrint opens files for printing parameter values at every step of execution. The file names are hard-coded, and all are stored within the subfolder files. Any existing files are over-written.
void SealCensus::StartPrint(void)
{
  hyperstream.setf(std::ios_base::scientific);
  hyperstream.precision(8);
  hyperfile="files/hyperParam.csv";    
  hyperstream.open(hyperfile.c_str());
  string header = "N \t Arrive \t SDArrive \t Correlation \t Tenure \t TenureCV ";
  hyperstream << "Step \t " << header << endl;

  hyperSDstream.setf(std::ios_base::scientific);
  hyperSDstream.precision(8);
  hyperfileSD="files/hyperSDParam.csv";    
  hyperSDstream.open(hyperfileSD.c_str());
  hyperSDstream << "Step \t " << header << "\t Like \n";

  for(auto it=allyear.begin();it!=allyear.end();it++)
   {
	string filenameYear="files/Year"+to_string(*it)+".csv";
	yearfile.insert(make_pair(*it,filenameYear));
	yearstream.open(filenameYear.c_str());
	yearstream << "Step \t " << header << "\t Like \n";
    yearstream.close();
   }
   
}

// Print one set of hyper parameters and hyperSD parameters. The streams hyperstream and hyperSDstream are held open throughout the run.
void SealCensus::PrintHyper()
{
	int j;
	hyperstream << counter << "\t";
	for(j=0; j<(Npar-1); j++) hyperstream << hyper[j] << "\t";
	hyperstream << hyper[Npar-1] << endl;
	hyperSDstream << counter << "\t";
	for(j=0; j<Npar; j++) hyperSDstream << hyperSD[j] << "\t";
	hyperSDstream << hyperLikelihood << endl;
}


// Print one set of parameters for all years. There is a file for each year, but a single file pointer is used for all of them. The pointer is opened to a single year's file, the parameters written to it, then the pointer closed. 
void SealCensus::PrintYears()
{
  for(auto it=allyear.begin();it!=allyear.end();it++)
   {
	yearstream.open(yearfile[*it].c_str(),ios_base::app);
	yearstream << counter << "\t";
	for(int j=0; j<Npar; j++) yearstream << fullparam[*it][j] << "\t";
	yearstream << fullLikelihood << endl;
    yearstream.close();
   }
}

// The class function ReturnSummary does the work of calculating summary statistics for the parameter chains stored in all the parameter files. It makes use of a class ReadParam (defined in readMetropParam.h). That class finds the mean and 95th percentiles (for lower and upper credible intervals) for each column in every one of the files. Those means and credible intervals are printed into the file ParameterResults.csv.
void SealCensus::ReturnSummary()
{
   int burnin=stepping[1];
   string summaryfile="ParameterResults.csv";
   output.setf(std::ios_base::scientific);
   output.precision(8);
   output.open(summaryfile.c_str());
   string header = "Year \t Stat \t N \t Arrive \t SDArrive \t Correlation \t Tenure \t TenureCV \n";
   output << header;
   
   for(auto it=allyear.begin();it!=allyear.end();it++)
    {
      ReadParam finalOneYear(yearfile[*it],burnin);
      unordered_map<string,NumVector> y=finalOneYear.CalcSummary();
      PrintOneResult(y,"Mean",to_string(*it));
      PrintOneResult(y,"Lower",to_string(*it));
      PrintOneResult(y,"Upper",to_string(*it));
    }

   ReadParam finalHyper(hyperfile,burnin);	
   ReadParam finalHyperSD(hyperfileSD,burnin);	
   unordered_map<string,NumVector> summaryHyper=finalHyper.CalcSummary();
   unordered_map<string,NumVector> summaryHyperSD=finalHyperSD.CalcSummary();

   PrintOneResult(summaryHyper,"Mean","Hyper");
   PrintOneResult(summaryHyper,"Lower","Hyper");
   PrintOneResult(summaryHyper,"Upper","Hyper");
   PrintOneResult(summaryHyperSD,"Mean","HyperSD");
   PrintOneResult(summaryHyperSD,"Lower","HyperSD");
   PrintOneResult(summaryHyperSD,"Upper","HyperSD");
   output.close();
}

// The class function PrintOneResult prints a single vector of summary parameters, either Mean, Lower, or Upper, for a single year, or for hyper-parameters.
void SealCensus::PrintOneResult(unordered_map<string,NumVector> result, string which, string yr)
{
  output << yr << "\t" << which << "\t";
  if(yr=="Hyper" || yr=="HyperSD") output << "\t";
  else output << result[which][1] << "\t";
  for(int i=2;i<6;i++) output << result[which][i] << "\t";
  output << result[which][6] << endl;
}

#endif
