/*
  Seal molt in Cpp
*/

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <iterator>
#include "/home/condit/programs/utilities/utilities.cpp"
#include "/home/condit/programs/utilities/likelihood.h"
#include "/home/condit/programs/utilities/statistics.h"
#include "/home/condit/programs/utilities/randomgenerator.h"
#include "/home/condit/elephantseal/phenology/cpp/modelMolt.h"

using namespace std;

/* Compilation and sample run. Folder name is optional.
cd ~/elephantseal/phenology/cpp
g++ -Wall modelMolt.cpp -o modelMolt.exe 2>&1 | more
g++ -Wall modelMolt.cpp -o modelMolt.exe

./modelMolt.exe /home/condit/elephantseal/phenology/ FemaleMolt 4000 500 1 test
*/

int main(int argc, char *argv[])
{
  string path=argv[1];
  string file=argv[2];
  int steps=atoi(argv[3]); 
  int show=atoi(argv[4]); 
  int test=atoi(argv[5]);
  string lab(argv[6]);
  
  MoltModel(path,file,steps,show,test,lab);
}


/* Options in R (not at github)
sealphenology(data=TRUE)
Testing smaller set:
prepareMoltCpp(moltdata=subset(MOLTDATA22,year>=2020))

label parameter for model runs, all 2016-2022--
30Nov2022-Run1: llikeYearAllHyper not included in updateGrandMean
30Nov2022-Run2: llikeYearAllHyper included in updateGrandMean
01Dec2022: llikeYearAllHyper included, fixed digits for sigma to 6
02Dec2022: core model, llikeYearAllHyper included, cleaned up molt data


*/


