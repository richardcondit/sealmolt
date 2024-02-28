# sealmolt
Software for fitting model to molt progress in elephant seals

Programs for fitting progress of molt in elephant seals, written in C++. THe program should work on a Linux system as is, with the data included, but knowledge of C++ and Bayesian statistics are to understand the output or adapt the code. The program fits a logistic curve to observations that increase from 0 to 100 for multiple individuals across multiple years. The logistic parameters include a hyperdistribution across years and animals.  

Files with extension .h or .cpp contain the code:
1) modelMolt.h defines the molt class and runs all calculations, and includes subroutines in header files
2) utilities.cpp has commonly used subroutines
3) util.h has commonly used subroutines
4) likelihood.h has common likelihood functions for statistical calculations
5) statistics.h has common statistical functions
6) randomgenerator.h has a random number generator
7) modelMolt.cpp has the main function, accepting command-line parameters and executing a complete run


The file FemaleMolt.csv has sample data. It includes 4 tab-delimited columns:
1) animalID is an arbitrary identifier of individuals, always an integer
2) year is an integer
3) yday is day of the year, with 1=1Jan, treated as a double
4) molt is observation of the molt status of the given animal on the given day, always in [0,100], a double

First 7 rows of data showing observed molt progress in animal 44168 in year 2016:
animalID	year	yday	molt
44168	2016	126	0
44168	2016	129	0
44168	2016	130	0
44168	2016	137	0
44168	2016	139	5
44168	2016	146	15
44168	2016	171	100


To illustrate compilation, assume file #7 and the data file are in a folder FOLDERNAME. The following creates an executable program in FOLDERNAME:
cd FOLDERNAME
g++ -Wall modelMolt.cpp -o modelMolt.exe 2>&1 | more
g++ -Wall modelMolt.cpp -o modelMolt.exe

Before executing, check lines 103-104 in modelMolt.h. They identify a folder where results will be written. Those should be renamed as needed, and both files must be created before executing. Then, to execute a run for 4000 steps, showing results every 500 steps, with output files named OUTPUT within those folders:
./modelMolt.exe FOLDERNAME FemaleMolt 4000 500 1 OUTPUT

There is one output file showing all parameter estimates of the Gibbs chain for every individual animal within the folder named pathFemale (line 104). There are files of parameters for each year and grand hyperparameters within the folder named path (line 103). Each file has 4000 rows for the 4000 steps chosen at run time.


