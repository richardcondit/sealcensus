#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H 1
 

#include <iostream>
#include <fstream>
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
// #include "/usr/lib/R/library/Rcpp/include/Rcpp.h"
#include <random>
#include "/home/condit/programs/utilities/util.h"


// This class provides easy functions for random draws to use from C++, based on the std random. The seeding, srand, must happen when initialized. 

class RandomGenerator
  {
    public:
      RandomGenerator();
      NumVector norm(int N, double mean, double sd);
      NumVector unif(int N, double low, double high);
      IntVector binom(int N, int S, double prob);
      IntVector pois(int N, double lambda);
      IntVector negbinom(int N, double mu, double k);
      NumVector gamma(int N, double shape, double scale);
  } ;

RandomGenerator::RandomGenerator() {  srand(time(NULL)); };
	
NumVector RandomGenerator::norm(int N, double mean, double sd)
{
 int i;
 vector<double> rvect;
 std::random_device rd;

 default_random_engine generator(rd());
 normal_distribution<double> distribution(mean,sd);

 for(i=0;i<N;i++) rvect.push_back(distribution(generator));
	
 return(rvect);
}

NumVector RandomGenerator::unif(int N, double low, double high)
{
 int i;
 vector<double> rvect;
 std::random_device rd;

 default_random_engine generator(rd());
 uniform_real_distribution<double> distribution(low,high);

 for(i=0;i<N;i++) rvect.push_back(distribution(generator));
	
 return(rvect);
}

IntVector RandomGenerator::binom(int N, int S, double prob)
{
 int i;
 IntVector rvect;
 std::random_device rd;

 default_random_engine generator(rd());
 binomial_distribution<int> distribution(S,prob);

 for(i=0;i<N;i++) rvect.push_back(distribution(generator));
	
 return(rvect);
}

IntVector RandomGenerator::pois(int N, double lambda)
{
 int i;
 IntVector rvect;
 std::random_device rd;

 default_random_engine generator(rd());
 poisson_distribution<int> distribution(lambda);

 for(i=0;i<N;i++) rvect.push_back(distribution(generator));
	
 return(rvect);
}

NumVector RandomGenerator::gamma(int N, double shape, double scale)
{
 int i;
 NumVector rvect;
 std::random_device rd;

 default_random_engine generator(rd());
 gamma_distribution<double> distribution(shape,scale);

 for(i=0;i<N;i++) rvect.push_back(distribution(generator));
	
 return(rvect);
}

// This fails if the first parameter of the distribution, which is size in R's rnbinom, is < 1. It does not fail in R. This does allow a non-integer size parameter, though, so if the first parameter is 1.001 it works. Unfortunately, very small values of size (ie k for clumping) are required for highly clumped distributions. 
IntVector RandomGenerator::negbinom(int N, double mu, double k)
{
 int i;
 IntVector rvect;
 double prob;
 
 prob=k/(k+mu);
 
 std::random_device rd;

 default_random_engine generator(rd());
 negative_binomial_distribution<int> distribution(k,prob);

 for(i=0;i<N;i++) rvect.push_back(distribution(generator));
	
 return(rvect);
}

#endif
