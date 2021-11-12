# sealcensus
SOFTWARE FOR ESTIMATING POPULATION SIZE OF A COLONY IN WHICH INDIVIDUALS ARE ASYNCHRONOUS
12 November 2021


This project archive includes all the files needed to compile and execute a program to analyze many daily counts of seals, or individuals of any other organism, across many seasons. The goal is estimating the total number of animals at a location in cases where animals are not all present at once. Details about the model and its application in elephant seals are presented in 

Condit, R., Le Boeuf, B. J., Morris, P. A., and Sylvan, M. 2007. Estimating population size in asynchronous aggregations: a Bayesian approach and test with elephant seal censuses. Marine Mammal Science 23:834–855. http://conditdatacenter.org/pdfs/ConditEtAlMMS2007.pdf.

Condit, R., Allen, S.G., Costa, D.P., Codde, S., Goley, P.D., Le Boeuf, B.J., Lowry, M.S., Morris, P. 2020. Estimating the breeding population of an elephant seal colony from a single census. bioRxiv 2020.12.01.368605; https://doi.org/10.1101/2020.12.01.368605.

Knowledge of hierarchical modeling and Bayesian estimation methods will be needed to understand this program, ie to change it. But it is designed to compile and execute without any such knowledge, given the following instructions.

THE ARCHIVE

There are 10 files in the archive, including this README.txt file, plus a folder named 'files'. Those with extension .h or .cpp are the program files. Those with .csv or .txt carry the input data needed to run the program. 

README.txt
sealcensusfull.cpp: main function for executing the full program
sealcensusFullCpp.h: definition of the SealCensus class, in which the model of population size is completed
readMetropParam.h: definition of the ReadParam class, in which tables Metropolis Markov-Chain parameter chains are converted to model output
Array2D.h: a utility defining a 2-dimensional data array (a cpp vector of vectors)
randomgenerator.h: a random number generator, needed for the Metropolis updates
statistics.h: many standard statistics functions, including Sum and pdfNorm used in SealCensus
utilities.cpp: many basic utility functions, including Test and TokenizeStr used in SealCensus
SampleData.csv: the file with all input data; this is a sample which can be used for testing the model, then to be replaced by user
inputfile.txt: a file with input parameters for the model; this is a sample relevant for elephant seals, to be adjusted as needed by user


INPUT COUNT DATA

Input data (as in SampleData.csv, but file can be renamed) must have 3 columns of tab-delimited integers

     -- There should be no header row
     
     -- Column 1 must be Season: an integer, in typical use, this is a year, but other numbers will work
     
     -- Column 2 must be Day within Season: an integer, in typical use, this is the day within each year
     
     -- Column 3 is the Count on each day: must be an integer
     
     -- There might be only one season, but then hyper-parameters are meaningless
     
     -- Days must have the same meaning each year, ie day 10 might mean 10 Jan every year (but can mean any other day)
     
     -- It is not necessary to have matching days every year; one year can have day 1, 5, 10, the next year 7, 8, 12
     
     -- Only include days with counts, no blank records; some seasons may have no data at all
     
     -- There can be seasons with few, or even one day, but some seasons need a full series of counts

INPUT MODEL PARAMETERS

Only 6 input parameters are needed from a file, and must be placed on two rows (as in inputfile.txt, but the file can be renamed). The first row has Start Day and End Day, the second row has 4 prior distributions.

  Start Day -- Start day each season should be earlier than counts ever start, ie before the earliest arrivals. It can be negative. The results are most accurate if there are days with predicted count~0. The default -30 works for northern elephant seals; it is equivalent to 31 October. Day 1 is then 1 December.

  End Day -- End day each season should be a day later than the last departure.

  Mean Tenure -- The mean length of time an individual is present. This must be known independently for good estimates. The default (in days) is for northern elephant seal.

  Prior SD of Tenure -- The prior standard deviation (SD) of that mean tenure. It is the degree of confidence in the independent estimate. Ideally, it is very small; if it is high, it will add error to the estimated population size. It must be positive. The default is from northern elephant seals.

  CV Tenure -- The coefficient of variation (CV) in tenure among individuals (CV = ratio SD/mean). This must be known independently for good estimates. Note the difference between CVtenure, which is a trait of the organism, and prior SD of tenure, which is confidence of the observations.

  Prior SD of CV Tenure -- The prior standard deviation (SD) of that CV tenure. It is the degree of confidence in the independent estimate of CV tenure. Ideally, it is very small; if it is high, it will add error to the estimated population size. It must be positive. The default is from northern elephant seals.


INSTRUCTIONS FOR COMPILATION

In Linux, this command from within the folder where all files were downloaded should successfully compile the executable:
  
  > g++ -Wall sealcensusFull.cpp -o sealcensus.exe

The randomgenerator will fail in old versions of gcc. The name sealcensus.exe can of course be changed.


INSTRUCTIONS FOR EXECUTION
  
Before executing, a subfolder named 'files' must be created within that main folder (containing the compiled sealcensus.exe). Its name is hard-coded in the program, so it must be present with that name. Once 'files' is created, from the same Linux folder, sealcensus.exe (or whatever the executable was named) must be executed with 5 command-line parameters. This one execution completes the entire analysis and prints results to a file. For example:
  
  > ./sealcensus.exe SampleData.csv inputfile.txt 6000 1000 999

  -- Argument 1 is the name of the data file, here the included SampleData.csv. A new data file should be created with a different name, then SampleData.csv can be preserved in case it is needed again for testing.
  
  -- Argument 2 is the name of the file with input parameters, here inputfile.txt. Again, a new file with a different name should be created for the particular case needed. See the descriptions of the 6 parameters needed given above.
  
  -- Argument 3 is the number of steps the model will run the parameter search. Final results should be based on 6000-10000 steps. Initial runs should be ~200 steps to confirm the model runs. A few hundred will finish quickly even with large datasets having many seasons. A final run of 6000 or more may take 10 minutes with many seasons.
  
  -- Argument 4 is the number of preliminary steps to be discarded in parameter calculations. Must be less than the number of steps. Final run should be 2000. For testing purposes, it matters little, as long as it is < steps.

  -- Argument 5 controls display to the screen while executing, eg if argument 5 is 999, the current estimated hyper-parameters will print to the screen every 999 steps. This helps confirm the model is working and estimates are converging. For testing, it might be 10 or 50, but for the final run 100 or more; if > steps, nothing will display until completion.
  
As the program runs, a completed set of parameters is printed to ascii files inside the folder named 'files'. There is one file for every season submitted, plus a set of hyper-parameters and their standard deviations. 

The results are printed to a file ParameterResults.csv in the same folder as the execution. The papers cited at the top should be consulted to understand details of the output. The output file name is hard-coded, so the user will need to copy the results file elsewhere to re-run the program with new data. The complete set of parameters, within the folder named 'files', should also be stored if needed by the user, because they will be over-written with any run. 

The output file, ParameterResults.csv, includes the estimated population parameters for every season submitted, plus the hyper-parameters across seasons. 

Users are welcome to contact me with questions. My email is in the publications listed at the top.


Richard Condit

