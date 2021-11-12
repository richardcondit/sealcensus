#ifndef STATISTICS_H
#define STATISTICS_H 1

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <iterator>
#include <math.h>
#include <map>
#include <set>
#include <algorithm>
#include "/home/condit/programs/utilities/utilities.cpp" 

// Basic math and probability functions, or generic string functions (StrSplit)

double Mean(NumVector x)
{
  double sum=0;
  int counter=0;
  for(size_t i=0;i<x.size();i++) if(!std::isnan(x[i]))
	 {
		 sum+=x[i];
		 counter++;
	 }
  if(counter>0) return(sum/x.size());
  return(nan(""));
}


double Mean(IntVector x)
{
	int counter=0;
    double sum=0;
	for(size_t i=0;i<x.size();i++) if(!std::isnan(x[i]))
	 {
		 sum+=x[i];
		 counter++;
	 }
	if(counter>0) return(sum/x.size());
	return(nan(""));
}


double Var(NumVector x)
{
	if(x.size()<2) return(0);
	double sum=0, mn=Mean(x);
	for(size_t i=0;i<x.size();i++) sum+=pow(x[i]-mn,2);
	return(sum/(x.size()-1));
}

double Var(IntVector x)
{
	if(x.size()<2) return(0);
    double sum, mn=Mean(x);
	for(size_t i=0;i<x.size();i++) sum+=pow(x[i]-mn,2);
	return(sum/(x.size()-1));
}

double SD(NumVector x)
{
	return(sqrt(Var(x)));
}

double SD(IntVector x)
{
	return(sqrt(Var(x)));
}

// Sort a vector to return a single quantile.
double Quantile(NumVector x, double prob)
{
  if(prob<=0) return(nan(""));
  if(prob>=1) return(nan(""));
  int len=x.size(), position;
  if(len==0) return(nan(""));
  
  NumVector sorted=x;
  sort(sorted.begin(),sorted.end());
  position=rint(len*prob);

  return(sorted[position]);
}

// Combinatorial based on recursion
int Choose(int N, int r)
{
 if(r<0 || r>N) return(0);
 if(r==1) return(N);
 if(r==N) return(1);	
 if(2*r<N) r=N-r;
 
 // Instead of recursion, perhaps it is faster to use logChoose, based on lgamma, then invert the log. But I seldom use Choose, only logChoose.
 return((int)(Choose(N-1,r-1)*N/r));
}

// Separate function for log of choose in case it is safer or faster for large numbers
double logChoose(int N, int r)
{
 if(r<0 || r>N) return(-INFINITY);
 if(r==1) return(log(N));  
 if(r==N) return(0);	

 // C++ has lgamma, which may or may not be faster than the recursive formula
 return(lgamma(N+1)-lgamma(r+1)-lgamma(N-r+1));
}


// Gaussian probability for a single x, logged or not based on boolean loglike.
double pdfNorm(double x, double mean, double sd, bool loglike)
{
	double denom=sd*sqrt(2*M_PI);
	double expon=pow((x-mean),2)/(2*sd*sd);
	
	if(!loglike) return(exp(-expon)/denom);
	
	return(-expon-log(denom));
}

// Gaussian probability density for a vector x, given constant mean and sd. Returns the full vector of probabilities, either logged or not. 
NumVector pdfNorm(NumVector x, double mean, double sd, bool loglike)
{
	NumVector prob;
	size_t i, len=x.size();
	
	for(i=0;i<len;i++) prob.push_back(pdfNorm(x[i],mean,sd,loglike));
	 
    return(prob);
}

// Version 3: Gaussian probability density for a vector x, with vector for mean and sd. Returns the full vector of probabilities, either logged or not. 
NumVector pdfNorm(NumVector x, NumVector mean, NumVector sd, bool loglike)
{
	NumVector prob;
	size_t i, len=x.size();
	NumVector bad(len,NAN);
	if(mean.size()!=len || sd.size()!=len) return(bad);
	
	for(i=0;i<len;i++) prob.push_back(pdfNorm(x[i],mean[i],sd[i],loglike));
	 
    return(prob);
}

// Binomial probability for a single atomic sample, logged or not based on boolean loglike, but calculated using logs.
double pdfBinom(int success, int trial, double prob, bool loglike)
{
	bool bail=0;
	double logC, logsucceed, logfail, logResult;
	
	if(std::isnan(success) || std::isnan(trial) || std::isnan(prob)) return(NAN);
	if(prob>1 || prob<0) return(NAN);
	if(success>trial || success<0) return(NAN);

	if(prob==1 && success==trial) { logResult=0; bail=1; }
	if(prob==1 && success<trial) { logResult=(-1)*INFINITY; bail=1; }
	if(prob==0 && success>0) { logResult=(-1)*INFINITY; bail=1; }
	if(prob==0 && success==0) { logResult=0; bail=1; }
	if(bail && loglike) return(logResult);
	if(bail && !loglike) return(exp(logResult));
	
	logC=logChoose(trial,success);
	logsucceed=success*log(prob);
	logfail=(trial-success)*log(1-prob);
    logResult=logC+logsucceed+logfail;
	
	if(!loglike) return(exp(logResult));
	
	return(logResult);
}

// Binomial probability density for a vector of observations, with atomic T and p. Returns the full vector of probabilities, either logged or not. 
NumVector pdfBinom(IntVector S, int T, double p, bool loglike)
{
	NumVector prob;
	size_t i, len=S.size();
	
	for(i=0;i<len;i++) prob.push_back(pdfBinom(S[i],T,p,loglike));
	 
    return(prob);
}

// Version 3: Binomial probability density for a vector of observations with matching vectors of T and p. Returns the full vector of probabilities, either logged or not. 
NumVector pdfBinom(IntVector S, IntVector T, NumVector p, bool loglike)
{
	NumVector prob;
	size_t i, len=S.size();
	NumVector bad(len,NAN);
	
	if(T.size()!=len || p.size()!=len) return(bad);
	
	for(i=0;i<len;i++) prob.push_back(pdfBinom(S[i],T[i],p[i],loglike));
	 
    return(prob);
}

// Poisson probability for a single atomic sample, logged or not based on boolean loglike, but calculated using logs.
double pdfPoisson(int obs, double lambda, bool loglike)
{
	double xFact, logNumer, logResult;
	
	if(std::isnan(obs) || std::isnan(lambda)) return(NAN);
	if(lambda<=0) return(NAN);
	if(obs<0) return(NAN);

	xFact=lgamma(obs+1);
	logNumer=obs*log(lambda)-lambda;
	logResult=logNumer-xFact;
	
	if(!loglike) return(exp(logResult));	
	return(logResult);
}

// Poisson probability density for a vector of observations, with atomic p. Returns the full vector of probabilities, either logged or not. 
NumVector pdfPoisson(IntVector S, double p, bool loglike)
{
	NumVector prob;
	int i, len=S.size();
	for(i=0;i<len;i++) prob.push_back(pdfPoisson(S[i],p,loglike));
    return(prob);
}

// Poisson probability density for a vector of observations with matching vectors of S and p. Returns the full vector of probabilities, either logged or not. 
NumVector pdfPoisson(IntVector S, NumVector p, bool loglike)
{
	NumVector prob;
	size_t i, len=S.size();
	NumVector bad(len,NAN);
	if(p.size()!=len) return(bad);
	for(i=0;i<len;i++) prob.push_back(pdfPoisson(S[i],p[i],loglike));
    return(prob);
}

// Exponential distribution for atomic observation and rate constant. Returns one probability, either logged or not. 
double pdfExp(double obs, double k, bool loglike)
{
	if(std::isnan(obs) || std::isnan(k)) return(NAN);
	if(k<=0) return(NAN);
	if(k<0) return(NAN);
	double logProb=(-1)*k*obs;
	if(!loglike) return(exp(logProb));	
	return(logProb);
}

// Exponential distribution for vector of observations and a single rate constant. Return the full vector of probabilities, either logged or not. 
NumVector pdfExp(NumVector S, double k, bool loglike)
{
	NumVector prob;
	size_t i, len=S.size();
	for(i=0;i<len;i++) prob.push_back(pdfExp(S[i],k,loglike));
    return(prob);
}

// Exponential distribution for vector of observations and matching rate constants. Returns the full vector of probabilities, either logged or not. 
NumVector pdfExp(NumVector S, NumVector k, bool loglike)
{
	NumVector prob;
	size_t i, len=S.size();
	NumVector bad(len,NAN);
	if(k.size()!=len) return(bad);
	for(i=0;i<len;i++) prob.push_back(pdfExp(S[i],k[i],loglike));
    return(prob);
}


// Find index of minimum. Must return as NumVector, because ties mean 2 or more are returned.
IntVector MinIndex(NumVector x)
{
	double lowest=x[0];
    IntVector indices;
	for(size_t i=0;i<x.size();i++) 
	   {
	     if(x[i]<lowest) { indices.clear(); lowest=x[i]; }
	     if(x[i]<=lowest) { indices.push_back(i); lowest=x[i]; }
	   } 
    return indices;
}

// Find index of maximum. Must return as NumVector, because ties mean 2 or more are returned.
IntVector MaxIndex(NumVector x)
{
	double highest=x[0];
    IntVector indices;
	for(size_t i=0;i<x.size();i++) 
	   {
	     if(x[i]>highest) { indices.clear(); highest=x[i]; }
	     if(x[i]>=highest) { indices.push_back(i); highest=x[i]; }
	   } 
    return indices;
}

// Copied from statisticsRcpp.h, removing IntVector, NumVector 
int Sum(vector<int> x)
{
	int total=0;
	for(auto it=x.begin();it<x.end();it++) total+=(*it);
	return(total);
}

double Sum(vector<double> x)
{
	double total=0;
	for(auto it=x.begin();it<x.end();it++) total+=(*it);
	return(total);
}

double Min(vector<double> x)
{
	double lowest=x[0];
	for(auto it=x.begin();it!=x.end();it++) if((*it)<lowest) lowest=(*it);
	return(lowest);
}

int Min(vector<int> x)
{
	int lowest=x[0];
	for(auto it=x.begin();it!=x.end();it++) if((*it)<lowest) lowest=(*it);
	return(lowest);
}

double Max(vector<double> x)
{
	double highest=x[0];
	for(auto it=x.begin();it!=x.end();it++) if((*it)>highest) highest=(*it);
	return(highest);
}

int Max(vector<int> x)
{
	int highest=x[0];
	for(auto it=x.begin();it!=x.end();it++) if((*it)>highest) highest=(*it);
	return(highest);
}

int Min(set<int> x)
{
	int lowest=(*x.begin());
	for(auto it=x.begin();it!=x.end();it++) if((*it)<lowest) lowest=(*it);
	return(lowest);
}

int Max(set<int> x)
{
	int highest=(*x.begin());
	for(auto it=x.begin();it!=x.end();it++) if((*it)>highest) highest=(*it);
	return(highest);
}


#endif
