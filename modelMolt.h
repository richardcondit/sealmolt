#ifndef MODELMOLT_CPP_H
#define MODELMOLT_CPP_H 1

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
#include "/home/condit/programs/utilities/utilities.cpp"
#include "/home/condit/programs/utilities/util.h"
#include "/home/condit/programs/utilities/likelihood.h"
#include "/home/condit/programs/utilities/statistics.h"
#include "/home/condit/programs/utilities/randomgenerator.h"
// #include "/home/condit/elephantseal/phenology/cpp/modelPart.h"

using namespace std;

/*
Moving the model of molt dates into cpp, following the birth model. This is a simplified model, Nov 2022, following two years of experience running an R version. This omits the non-detection probabilities, and omits the linear change through time. See modelMoltOrig.h with full original version, but likely obsolete.

Data format match birth model except of course the column pctmolt instead of pup:
* A set female has a vector animalIDs, iterated to work through all females
* A map femaleyr with animalID as key has a set of unique years for each female, iterated to work through one female's years 
* A 2-level map X has observations. The outer level is for individual females, with animalID as key. For each, there is an inner map of years, with key=year and each element a Array2D of double, 2 columns, year, yday, molt. So X[1909][2022] accesses data as Array2D for female 1909 in year 2022. 

Parameter format are very close to birth model, with 2-parameter logistic replacing the 2-parameter birth-loss function. The other conspicuous difference is that the molt model requires a SD for the logistic, and it will be the same across all females and years -- a fixed effect. Also distinct from the birth model, I am omitting a regression between molt day and time. 
* The hyper-distribution for molt day is the complex one. My goal is both a year term and a female term, requiring many SDs and scale parameters.
* The hyper-distribution for logisitc slope will differ only in that I will ignore annual terms.
* grand mean of all molt days, across females and years 
  -- grandmean double
* map for mean molt day each year 
  -- yearmean Int2Num, with year as key and each element a double
* map for lifetime mean molt day for each female 
  -- femalemean Int2Num, with animalID as key and each element a double
* map for lifetime mean slope (of logistic) for each female 
  -- femalemeanslope Int2Num, with animalID as key and each element a double
* map for vector of annual molt day per female
  -- femalemolt Int2Vect, a map of maps, animalID as outer key, year as inner, a molt day for each year
* map for vector of annual logistic slope per female
  -- femaleslope Int2Vect, a map of maps, animalID as outer key, year as inner, a molt day for each year
* sigma for error of molt function
  -- sigma double
* SD within year, multiple molt days years within a year
  -- double withinyearSD, identical for all years
* SD within female, multiple years within a female
  -- double withinfemaleSD, identical for all females
* SD among yearmean, multiple annual means
  -- double amongyearSD
* SD among femalemean, multiple female lifetime
  -- double amongfemaleSD

* Scale parameters match exactly
* 
*/

struct MoltRecord { int animalID; int year; double yday; double molt; };
typedef vector<MoltRecord> MoltRecords;
typedef map<int,Array2D<NumVector>> MoltGroup;
typedef map<int,MoltGroup> MoltBiGroup;
typedef map<int,set<int>> Int2Int;
typedef map<int,double> Int2Num;
typedef map<int,Int2Num> Int2Vect;

class MoltModel
{    
  protected:
    RandomGenerator rand;
    int counter=0, N=0, test;
    bool logit=TRUE;
    double target=0.25, adjexp, adjust=1.01, startscale;
    IntVector stepping;

    // Observations
    MoltRecords x;               // Complete data table, no grouping
    MoltBiGroup xByFemaleYr;     // Complete data, grouped at two levels: outer animalID, inner year
    set<int> female;             // AnimalIDs of all females
    set<int> year;               // All years
    Int2Int femaleyr;            // All years for each female 

    // Parameters
    double sigma=0.1, grandmean=130.0, withinyearSD=10.0, withinfemaleSD=1.0, amongyearSD=2.0, amongfemaleSD=10.0,
           grandmeanslope=6.0, withinfemaleslopeSD=1.0, amongfemaleslopeSD=1.0;
    Int2Num yearmean, femalemean, femalemeanslope;
    Int2Vect femalemolt, femaleslope;

    // Scale
    double scaleslope=0.6, scalesigma=0.01, scalegrand=5.0, scalewithinyearSD=1.0, scalewithinfemaleSD=0.1, scaleamongyearSD=0.2, scaleamongfemaleSD=1.0,
           scalegrandmeanslope=0.6, scalewithinfemaleslopeSD=0.1, scaleamongfemaleslopeSD=0.1;
    Int2Num scaleyearmean, scalefemalemean, scalefemalemeanslope;
    Int2Vect scalefemalemolt, scalefemaleslope;
  
    // Full likelihood each step
    double stepLikelihood;
    
    // Output files
    string path="/home/condit/elephantseal/phenology/cpp/lmerMolt/";
    string pathFemale="/home/condit/elephantseal/phenology/cpp/lmerMoltFemale/";
    FILE *mainfile, *moltfile, *yearfile;
    string testID, outputname, extension=".csv", suffix;

  public:
    MoltModel(string path, string infile, int steps, int show, int debug, string label);
    void Initialize(string);
    void FillMaps(void);
    void StartParam(void);
    double llikeMolt(int,int,double,double,double);
    double llikeAllMolt(double,bool);
    double llikeFemaleLifetime(Int2Num,double,double,double,double);
    double llikeMoltWithinYear(int,double,double,double,double);
    double llikeLifetimeAllHyper(double,double,int);
    double llikeYearAllHyper(double,double);
    double llikeLifetimeObs(Int2Num,double,double);
    double llikeLifetimeAllObs(double);
    double llikeLifetimeAllObsSlope(double);
    double llikeYearAllObs(double);
    double updateOneMid(int,int);
    double updateOneSlope(int,int);
    double updateSigma(double);
    double updateFemaleMean(int);
    double updateFemaleMeanSlope(int);
    double updateYearMean(int);
    double updateGrandMean(double);
    double updateAmongSD(double);
    double updateAmongYearSD(double);
    double updateWithinSD(double);
    double updateWithinYearSD(double);
    double updateGrandMeanSlope(double);
    double updateAmongFemaleSlopeSD(double);
    double updateWithinFemaleSlopeSD(double);
    void updateFemales(void);
    void updateYears(void);
    void FullRun(void);
    NumVector logisticMolt(NumVector, NumVector);
    void StartPrint(void);
    void PrintMain(void);
    void PrintFullFemale(void);
    void PrintFullYear(void);
};

// After path and data file name, the 4 integers are steps, show, and a debug flag
MoltModel::MoltModel(string path, string infile, int steps, int show, int debug, string label)
{
  test=debug;
  suffix=label;
  testID=to_string(test);
  stepping.push_back(steps);
  stepping.push_back(0);      // burn-in is not used within Cpp
  stepping.push_back(show);  
  adjexp=(1-target)/target;
    
  string folder="cpp/";
  string outputname=path+folder+extension;
  
  Initialize(path+folder+infile+extension);
  FullRun();
}

void MoltModel::Initialize(string filename)
{
  counter=0;
  if(test==3) cout << filename << endl;
  adjexp=(1-target)/target;

  ifstream infile(filename);
  string oneline;
  getline(infile,oneline);   // Header row

  MoltRecord datum;
  set<int> emptyYr;
  MoltGroup emptyAnimal;
    
  while(getline(infile,oneline))
   {
    StrVector Cells;    
    TokenizeStr(oneline, Cells, "\t");
    
    datum.animalID=stoi(Cells[0]);
    datum.year=stoi(Cells[1]);
    datum.yday=stod(Cells[2]);
    datum.molt=stod(Cells[3])/100;
    
    x.push_back(datum);
    xByFemaleYr.insert(make_pair(datum.animalID,emptyAnimal));
    female.insert(datum.animalID);
    year.insert(datum.year);
    femaleyr.insert(make_pair(datum.animalID,emptyYr));
    N++;
    if(test==2) printf("Record %d is animal %d in year=%d, day %5.0lf, molt %6.3lf\n", N, datum.animalID, datum.year, datum.yday, datum.molt);
   }
  infile.close();

  FillMaps();
  StartParam();
  if(test==1) printf("%d records were loaded\n", N);
  if(test==1) Test(year);
  if(test==3) cout << "Number of females = " << xByFemaleYr.size() << endl;
}


// Filling all data objects, as well as the intercept-slope parameter for each female. Borrowed from lmerBayes but not identical since females are integers. But the 3-layer map of maps for data is considerably more complicated.
void MoltModel::FillMaps()
{
   // NumVector emptyParam(startparam), startscale{1.0,1.0,1.0};
   Array2D<NumVector> emptyData(0,2);

   // First fill in all years for every female
   for(auto it=x.begin();it!=x.end();it++) 
    {
      int ID=(*it).animalID;
      femaleyr[ID].insert((*it).year);
    }
    
   // Then create an empty table of molt observations, one per female-year, in variable xByFemaleYr
   for(auto fit=female.begin();fit!=female.end();fit++)
    {
       int ID=(*fit);
       set<int> herYears=femaleyr[ID];
       for(auto it=herYears.begin();it!=herYears.end();it++) 
        {
          int oneyr=(*it);
          xByFemaleYr[ID].insert(make_pair(oneyr,emptyData));
        }
    }
    
   // Finally fill all those tables per female-year
   for(auto it=x.begin();it!=x.end();it++) 
    {
      int ID=(*it).animalID;
      int yr=(*it).year;
      NumVector datum{(*it).yday,(*it).molt};
      xByFemaleYr[ID][yr].AddOneRow(datum);
    }
}   

// Fill maps of parameters and scale parameters: Int2Num yearmean, femalemean; Int2Vect femalemolt;
void MoltModel::StartParam(void)
{
   for(auto fit=female.begin();fit!=female.end();fit++)
    {
       int ID=(*fit);
       femalemean[ID]=130.0;
       femalemeanslope[ID]=6.0;
       scalefemalemean[ID]=10.0;
       scalefemalemeanslope[ID]=0.6;
       set<int> herYears=femaleyr[ID];
       for(auto it=herYears.begin();it!=herYears.end();it++) 
        {
          int oneyr=(*it);
          femalemolt[ID].insert(make_pair(oneyr,130.0));
          femaleslope[ID].insert(make_pair(oneyr,6.0));
        }
    }
    
   for(auto yit=year.begin();yit!=year.end();yit++)
    {
       int yr=(*yit);
       yearmean[yr]=130.0;
       scaleyearmean[yr]=1.0;
    }
}

void MoltModel::FullRun(void)
{
    StartPrint();
    
    for(int i=0;i<stepping[0];i++)
     {
        updateFemales();
        updateSigma(sigma);
        updateWithinSD(withinfemaleSD);
        updateAmongSD(amongfemaleSD);
        updateWithinFemaleSlopeSD(withinfemaleslopeSD);
        updateAmongFemaleSlopeSD(amongfemaleslopeSD);
        
        updateYears();
        updateAmongYearSD(amongyearSD);
        updateWithinYearSD(withinyearSD);

        updateGrandMean(grandmean);               
        updateGrandMeanSlope(grandmeanslope);               
        stepLikelihood=llikeAllMolt(sigma,1);

        PrintMain();
        PrintFullFemale();
        PrintFullYear();
        
        if((test<5 || test>10) && counter%stepping[2]==0) 
          {
             cout << "At step " << i << ", Likelihood: " << stepLikelihood;
             time_t rightnow = time(NULL);
             printf(" at time %s", ctime(&rightnow));
          }
        counter++;
        if(test==5 && counter%stepping[2]==0)
          {
             time_t rightnow = time(NULL);
             printf("Step %4.0d: sigma=%5.4lf, llike=%6.1lf, ", counter, sigma, stepLikelihood);
             printf("grandmean=%5.1lf, amongSD=%4.2lf, withinSD=%4.2lf, ", grandmean, amongfemaleSD, withinfemaleSD);
             printf("grandslope=%4.2lf, amongslopeSD=%4.3lf, withinslopeSD=%4.3lf ", grandmeanslope, amongfemaleslopeSD, withinfemaleslopeSD);
             printf("at %s", ctime(&rightnow));
          } 
     }

    fclose(mainfile);
    fclose(yearfile);
}

// Update the grand mean. This requires the female means with amongfemale SD, plus year means and amongyearSD. Either female or year means alone should lead to the same answer, so the question arises whether including both makes sense. 
double MoltModel::updateGrandMean(double mu)
{
    double origlike=llikeLifetimeAllHyper(mu,amongfemaleSD,0)+llikeYearAllHyper(mu,amongyearSD);  // Year component of likelihood
    double testval=rand.norm(1,mu,scalegrand)[0];
    if(testval<=0) return((-1)*INFINITY);
    double newlike=llikeLifetimeAllHyper(testval,amongfemaleSD,0)+llikeYearAllHyper(testval,amongyearSD);

    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)               // Accept
     {
        scalegrand=pow(adjust,adjexp);
        grandmean=testval;
        return(newlike);
     }
    scalegrand=(1/adjust);           // Reject
    return(origlike);
}

// Update the grand mean of female slopes. 
double MoltModel::updateGrandMeanSlope(double mu)
{
    double origlike=llikeLifetimeAllHyper(mu,amongfemaleslopeSD,1);
    double testval=rand.norm(1,mu,scalegrandmeanslope)[0];
    if(testval<=0) return((-1)*INFINITY);
    double newlike=llikeLifetimeAllHyper(testval,amongfemaleslopeSD,1);

    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)               // Accept
     {
        scalegrandmeanslope=pow(adjust,adjexp);
        grandmeanslope=testval;
        return(newlike);
     }
    scalegrandmeanslope=(1/adjust);           // Reject
    return(origlike);
}

// Update SD among female means. This is relevant only within years, so depends on all female molts within each year.
double MoltModel::updateAmongSD(double sd)
{
    double origlike=llikeLifetimeAllHyper(grandmean,sd,0);
    double testval=rand.norm(1,sd,scaleamongfemaleSD)[0];
    if(testval<=0) return((-1)*INFINITY);
    double newlike=llikeLifetimeAllHyper(grandmean,testval,0);

    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                       // Accept
     {
        scaleamongfemaleSD=pow(adjust,adjexp);
        amongfemaleSD=testval;
        return(newlike);
     }
    scaleamongfemaleSD=(1/adjust);           // Reject
    return(origlike);
}

double MoltModel::updateAmongFemaleSlopeSD(double sd)
{
    double origlike=llikeLifetimeAllHyper(grandmeanslope,sd,1);
    double testval=rand.norm(1,sd,scaleamongfemaleslopeSD)[0];
    if(testval<=0) return((-1)*INFINITY);
    double newlike=llikeLifetimeAllHyper(grandmeanslope,testval,1);

    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                       // Accept
     {
        scaleamongfemaleslopeSD=pow(adjust,adjexp);
        amongfemaleslopeSD=testval;
        return(newlike);
     }
    scaleamongfemaleslopeSD=(1/adjust);           // Reject
    return(origlike);
}

// SD among observations of each female over her lifetime. This depends on lifetime molt days and a mean for every female.
double MoltModel::updateWithinSD(double sd)
{
    double origlike=llikeLifetimeAllObs(sd);
    double testval=rand.norm(1,sd,scalewithinfemaleSD)[0];
    if(testval<=0) return((-1)*INFINITY);
    double newlike=llikeLifetimeAllObs(testval);

    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                        // Accept
     {
        scalewithinfemaleSD=pow(adjust,adjexp);
        withinfemaleSD=testval;
        return(newlike);
     }
    scalewithinfemaleSD=(1/adjust);           // Reject
    return(origlike);
}

// SD among observations of each female slopes over her lifetime. This depends on lifetime slopes and a mean for every female.
double MoltModel::updateWithinFemaleSlopeSD(double sd)
{
    double origlike=llikeLifetimeAllObsSlope(sd);
    double testval=rand.norm(1,sd,scalewithinfemaleslopeSD)[0];
    if(testval<=0) return((-1)*INFINITY);
    double newlike=llikeLifetimeAllObsSlope(testval);

    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                        // Accept
     {
        scalewithinfemaleslopeSD=pow(adjust,adjexp);
        withinfemaleslopeSD=testval;
        return(newlike);
     }
    scalewithinfemaleslopeSD=(1/adjust);           // Reject
    return(origlike);
}

// Update SD among year means. 
double MoltModel::updateAmongYearSD(double sd)
{
    double origlike=llikeYearAllHyper(grandmean,sd);
    double testval=rand.norm(1,sd,scaleamongyearSD)[0];
    if(testval<=0) return((-1)*INFINITY);
    double newlike=llikeYearAllHyper(grandmean,testval);

    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                     // Accept
     {
        scaleamongyearSD=pow(adjust,adjexp);
        amongyearSD=testval;
        return(newlike);
     }
    scaleamongyearSD=(1/adjust);           // Reject
    return(origlike);
}

// SD among observations of molt within years. This depends on all molt days and a mean, combining all years.
double MoltModel::updateWithinYearSD(double sd)
{
    double origlike=llikeYearAllObs(sd);
    double testval=rand.norm(1,sd,scalewithinyearSD)[0];
    if(testval<=0) return((-1)*INFINITY);
    double newlike=llikeYearAllObs(testval);

    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                        // Accept
     {
        scalewithinyearSD=pow(adjust,adjexp);
        withinyearSD=testval;
        return(newlike);
     }
    scalewithinyearSD=(1/adjust);           // Reject
    return(origlike);
}


// Loop through all females, and all years within each, updating molt date for each one.
void MoltModel::updateFemales(void)
{
  for(auto ity=femaleyr.begin();ity!=femaleyr.end();ity++)  
    {
      int f=ity->first;
      set<int> herYrs=ity->second;
      for(auto itf=herYrs.begin();itf!=herYrs.end();itf++)
       {
          int yr=(*itf);
          updateOneMid(f,yr);
          updateOneSlope(f,yr);
       }
      updateFemaleMean(f); 
      updateFemaleMeanSlope(f); 
     
      if(test==f && counter%stepping[2]==0)
       {
         printf("Female %6.0d: ", f);
         printf("mean(molt)=%5.1lf, mean(slope)=%4.2lf, ", femalemean[f], femalemeanslope[f]);
         printf("her Molts:\n"); 
         Test(femalemolt[f]);
         printf("her Slopes:\n"); 
         Test(femaleslope[f]);
       }
   }
}

// Update femalemolt, the molt day for a single female in a single year. This requires hyper-llikelihood!
double MoltModel::updateOneMid(int f, int y)
{
    double current=femalemolt[f][y];
    double hyper=pdfNorm(current,femalemean[f],withinfemaleSD,logit);
    double origlike=hyper+llikeMolt(f,y,current,femaleslope[f][y],sigma);
    double testval=rand.norm(1,current,scalefemalemolt[f][y])[0];
    hyper=pdfNorm(testval,femalemean[f],withinfemaleSD,logit);
    double newlike=hyper+llikeMolt(f,y,testval,femaleslope[f][y],sigma);
     
    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                             // Accept
     {
        scalefemalemolt[f][y]=pow(adjust,adjexp);
        femalemolt[f][y]=testval;
        return(newlike);
     }
    scalefemalemolt[f][y]=(1/adjust);              // Reject
    return(origlike);
}

// Exact parallel to updateOneMid, this to update femaleslope, the slope for a single female in a single year. This requires hyper-llikelihood!
double MoltModel::updateOneSlope(int f, int y)
{
    double current=femaleslope[f][y];
    double hyper=pdfNorm(current,femalemeanslope[f],withinfemaleslopeSD,logit);
    double origlike=hyper+llikeMolt(f,y,femalemolt[f][y],current,sigma);
    double testval=rand.norm(1,current,scalefemaleslope[f][y])[0];
    hyper=pdfNorm(testval,femalemeanslope[f],withinfemaleslopeSD,logit);
    double newlike=hyper+llikeMolt(f,y,femalemolt[f][y],testval,sigma);
     
    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                             // Accept
     {
        scalefemaleslope[f][y]=pow(adjust,adjexp);
        femaleslope[f][y]=testval;
        return(newlike);
     }
    scalefemaleslope[f][y]=(1/adjust);              // Reject
    return(origlike);
}

// Loop through all years, updating molt date for each one.
void MoltModel::updateYears(void)
{
  for(auto ity=year.begin();ity!=year.end();ity++)  
    {
      int yr=(*ity);
      updateYearMean(yr); 

      if(test==yr && counter%stepping[2]==0)
         printf("Step %5d: Year %4d molt=%6.1f, grandmean=%6.1lf, amongSD=%6.3lf, withinSD=%6.3lf, sigma=%5.4lf\n", 
                 counter, yr, yearmean[yr], grandmean, amongyearSD, withinyearSD, sigma);
    }
}

// Update the annual mean for a single year, given all female molts that year, plus hyper-llikelihood.
double MoltModel::updateYearMean(int y)
{
    double current=yearmean[y];
    // double hyper=pdfNorm(current,grandmean,amongyearSD,logit); Fuck-up: hyper is embedded within llikeMoltWithinYear!!
    double origlike=llikeMoltWithinYear(y,current,withinyearSD,grandmean,amongyearSD);
    double testval=rand.norm(1,current,scaleyearmean[y])[0];
    // hyper=pdfNorm(testval,grandmean,amongyearSD,logit);
    double newlike=llikeMoltWithinYear(y,testval,withinyearSD,grandmean,amongyearSD);
     
    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                         // Accept
     {
        scaleyearmean[y]=pow(adjust,adjexp);
        yearmean[y]=testval;
        return(newlike);
     }
    scaleyearmean[y]=(1/adjust);              // Reject
    return(origlike);
}

// Update the global error, sigma. It is fixed across years and females, so depends on all molt observations. Does not depend on any hyper-parameters.
double MoltModel::updateSigma(double sd)
{
    double origlike=llikeAllMolt(sd,0);
    double testval=rand.norm(1,sd,scalesigma)[0];
    double newlike=llikeAllMolt(testval,0);

    if(test==9 && counter%stepping[2]==0) 
      printf("Step %8d: orig=%6.3lf with sigma %6.3lf, new=%6.3lf with %6.3lf\n", counter, origlike, sd, newlike, testval);

    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                // Accept
     {
        scalesigma*=pow(adjust,adjexp);
        sigma=testval;
        return(origlike);
     }
    scalesigma*=(1/adjust);           // Reject
    return(newlike);
}

// Update the lifetime mean molt date for a single female. This depends on the grandmean of all female-years and its hyper-sigma, amongfemaleSD, as well as her observed molt dates.
double MoltModel::updateFemaleMean(int f)
{
    double current=femalemean[f];
    double origlike=llikeFemaleLifetime(femalemolt[f],current,withinfemaleSD,grandmean,amongfemaleSD);
    double testval=rand.norm(1,current,scalefemalemean[f])[0];
    double newlike=llikeFemaleLifetime(femalemolt[f],testval,withinfemaleSD,grandmean,amongfemaleSD);
    
    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                           // Accept
     {
        scalefemalemean[f]*=pow(adjust,adjexp);
        femalemean[f]=testval;
        return(newlike);
     }
    scalefemalemean[f]*=(1/adjust);           // Reject
    return(origlike);
}

// Update the lifetime mean slope for a single female. This depends on the grandmeanslope of all female-years and its hyper-sigma, amongfemaleslopeSD, as well as her observed slopes each year.
double MoltModel::updateFemaleMeanSlope(int f)
{
    double current=femalemeanslope[f];
    double origlike=llikeFemaleLifetime(femaleslope[f],current,withinfemaleslopeSD,grandmeanslope,amongfemaleslopeSD);
    double testval=rand.norm(1,current,scalefemalemeanslope[f])[0];
    double newlike=llikeFemaleLifetime(femaleslope[f],testval,withinfemaleslopeSD,grandmeanslope,amongfemaleslopeSD);
    
    if(std::isnan(newlike) && !std::isnan(origlike)) newlike=(-1)*INFINITY;
    double likeratio=exp(newlike-origlike);
    NumVector r=rand.unif(1,0.0,1.0);    
    if(r[0]<likeratio)                           // Accept
     {
        scalefemalemeanslope[f]*=pow(adjust,adjexp);
        femalemeanslope[f]=testval;
        return(newlike);
     }
    scalefemalemeanslope[f]*=(1/adjust);           // Reject
    return(origlike);
}

// The likelihood for slope and sigma parameters must use every single female observation. This works through every female in every year and calculates likelihood. With full set to TRUE, it adds the likelihood of parameters given hyperparameters.
double MoltModel::llikeAllMolt(double sd, bool full)
{
  double moltlike=0, modellike=0;
  if(sd<=0) return((-1)*INFINITY);
  
  for(auto ity=femaleyr.begin();ity!=femaleyr.end();ity++)  
    {
      int female=ity->first;
      set<int> herYrs=ity->second;
      
      for(auto itf=herYrs.begin();itf!=herYrs.end();itf++)
       {
          int yr=(*itf);
          moltlike+=llikeMolt(female,yr,femalemolt[female][yr],femaleslope[female][yr],sd);
       }

      if(full) modellike+=llikeFemaleLifetime(femalemolt[female],femalemean[female],withinfemaleSD,grandmean,amongfemaleSD);   
      if(full) modellike+=llikeFemaleLifetime(femaleslope[female],femalemeanslope[female],withinfemaleslopeSD,grandmeanslope,amongfemaleslopeSD);   
   }   

  return(moltlike+modellike);
}

// Likelihood over a female's entire lifetime depends on the estimated molts (slopes) in each year, her lifetime mean, with withinfemale(slope)SD, and the grandmean(slope) with amongfemale(slope)SD. This is used for updating various hyper-parameters, so all are passed: aSD is amongfemaleSD, wSD is withinfemaleSD, mu is grandmean, fmean is femalemean. 
double MoltModel::llikeFemaleLifetime(Int2Num allParam, double fmean, double wSD, double mu, double aSD)
{
    double hyperlike=pdfNorm(fmean,mu,aSD,logit);
    // Int2Num allMolts=femalemolt[f];
    double obslike=llikeLifetimeObs(allParam,fmean,wSD);
    return(obslike+hyperlike);
}

// Likelihood of all molts in one year depends on the estimated molts in that year, the year mean, with withinyearSD, and the grandmean with amongyearSD. This is used for updating various hyper-parameters, so all are passed: aSD is amongyearSD, wSD is withinyearSD, mu is grandmean, ymean is yearmean. 
double MoltModel::llikeMoltWithinYear(int y, double ymean, double wSD, double mu, double aSD)
{
    double hyperlike=pdfNorm(ymean,mu,aSD,logit);
    Int2Num allMolts;
    for(auto it=femalemolt.begin();it!=femalemolt.end();it++) 
     {
       int f=it->first;
       Int2Num onefemale=it->second;
       if(onefemale.find(y)!=onefemale.end()) allMolts.insert(make_pair(f,onefemale[y]));
       if(test==f && counter%stepping[2]==0)
           { cout << "At step " << counter << " Female " << f << " has molts: " << endl; Test(onefemale); }
     }
    if(test==9 && counter%stepping[2]==0)
           { cout << "At step " << counter << " year " << y << " has molts: " << endl; Test(allMolts); }
    double obslike=llikeLifetimeObs(allMolts,ymean,wSD);
    return(obslike+hyperlike);
}

// Loop through all females and get likelihood of their lifetime means, given amongfemaleSD and the grandmean. With which=0, its molt day, which=1 is slope
double MoltModel::llikeLifetimeAllHyper(double mu, double aSD, int which)
{
    double total=0;
    double herMean;
    for(auto it=female.begin();it!=female.end();it++) 
     {
        int ID=(*it);
        if(which==0) herMean=femalemean[ID];
        else herMean=femalemeanslope[ID];
        total+=pdfNorm(herMean,mu,aSD,logit);
     }
   return(total);
}

// Loop through all years and get likelihood of their means, given amongyearSD and the grandmean. There is no year variation in the slope, so no parameter for slopes in each yearS.
double MoltModel::llikeYearAllHyper(double mu, double aSD)
{
    double total=0;
    for(auto it=year.begin();it!=year.end();it++) 
     {
        int yr=(*it);
        double yrMean=yearmean[yr];
        total+=pdfNorm(yrMean,mu,aSD,logit);
     }
   return(total);
}


// Likelihood of a female annual molt dates (or slopes) relative to her own mean and withinfemaleSD (withinfemaleslopeSD).
double MoltModel::llikeLifetimeObs(Int2Num annual, double fmean, double wSD)
{
    NumVector obs;
    for(auto it=annual.begin();it!=annual.end();it++) obs.push_back(it->second);
    double obslike=Sum(pdfNorm(obs,fmean,wSD,logit));
    return(obslike); 
}

// Loop through all females and get likelihood of their annual molt dates relative to their lifetime means, given a withinfemaleSD.
double MoltModel::llikeLifetimeAllObs(double wSD)
{
    double total=0;
    for(auto it=female.begin();it!=female.end();it++) 
     {
        int ID=(*it);
        double herMean=femalemean[ID];
        Int2Num herMolts=femalemolt[ID];
        total+=llikeLifetimeObs(herMolts,herMean,wSD);
     }
   return(total);
}

// Loop through all females and get likelihood of their annual slopes relative to their lifetime means, given a withinfemaleslopeSD.
double MoltModel::llikeLifetimeAllObsSlope(double wSD)
{
    double total=0;
    for(auto it=female.begin();it!=female.end();it++) 
     {
        int ID=(*it);
        double herMean=femalemeanslope[ID];
        Int2Num herSlopes=femaleslope[ID];
        total+=llikeLifetimeObs(herSlopes,herMean,wSD);
     }
   return(total);
}

// Loop through all years and get likelihood of its molt dates relative to the year means, given a withinyearSD.
double MoltModel::llikeYearAllObs(double wSD)
{
    double total=0;
    for(auto it=year.begin();it!=year.end();it++) 
     {
        int yr=(*it);
        double yrMean=yearmean[yr];
        Int2Num yearMolts;
        for(auto itf=femalemolt.begin();itf!=femalemolt.end();itf++) 
         {
           int ID=itf->first;
           Int2Num onefemale=itf->second;
           if(onefemale.find(yr)!=onefemale.end()) yearMolts.insert(make_pair(ID,onefemale[yr]));
         }
        total+=llikeLifetimeObs(yearMolts,yrMean,wSD);
        // if(test==yr && counter%stepping[2]==0) Test(yearMolts);
     }
   return(total);
}

// Likelihood function for logistic model given observations of one female in one year and the two parameters for logistic. Needs a female and year to work with one at a time, but all parameters are passed so they can be updated.
double MoltModel::llikeMolt(int f, int y, double m, double s, double sd)
{
  double badlike=(-1)*INFINITY;  
  if(s<=0.01) return(badlike);
  if(m<60) return(badlike);

  Array2D<NumVector> obs=xByFemaleYr[f][y];  
  NumVector molt=obs.GetOneCol(1);
  NumVector day=obs.GetOneCol(0);
  NumVector par{m,s};
  NumVector sdvect(day.size(),sd);

  NumVector predmolt=logisticMolt(day,par);

  if(test<0 && y==2021 && counter%stepping[2]==0)
   {
     cout << "Female " << f << ": " << endl;
     Test(par);
     Test(sdvect);
     Test(predmolt);
   }
  
      
  bool logit=TRUE;
  NumVector llike=pdfNorm(molt,predmolt,sdvect,logit);

  return(Sum(llike));
}

// Open the files to print parameters. Each female needs two files, one for her model parameters and one for individual birth dates. One file covers hyperparameters and detection parameters. 
void MoltModel::StartPrint(void)
{
  string filename=path+"mainParam-"+suffix+".csv";    
  mainfile = fopen(filename.c_str(),"w");
  fprintf(mainfile,"counter \t amongfemaleSD \t withinfemaleSD \t amongfemaleslopeSD \t withinfemaleslopeSD \t ");
  fprintf(mainfile,"amongyearSD \t withinyearSD \t sigma \t grandmean \t grandmeanslope \t stepLikelihood \n");

  string filenameYear=path+"Year-"+suffix+".csv";
  yearfile=fopen(filenameYear.c_str(),"w");
  fprintf(yearfile,"step \t year \t molt \n");

  for(auto it=female.begin();it!=female.end();it++)
   {
    string ID=to_string(*it);
    string filenameFemale=pathFemale+suffix+"-animal-"+ID+".csv";
    moltfile=fopen(filenameFemale.c_str(),"w");
    fprintf(moltfile,"step \t animalID \t year \t molt \t meanmolt \t slope \t meanslope \n");
    fclose(moltfile);
   }
}

void MoltModel::PrintMain(void)
{
    fprintf(mainfile,"%6d \t %6.3lf \t %6.3lf \t %6.3lf \t %6.3lf \t %6.3lf \t %6.3lf \t %8.6lf \t %6.3lf \t %6.3lf \t %6.3lf \n", 
             counter, amongfemaleSD, withinfemaleSD, amongfemaleslopeSD, withinfemaleslopeSD, 
                      amongyearSD, withinyearSD, sigma, grandmean, grandmeanslope, stepLikelihood);
}

// Called at each step, print every female's current molt estimates.
void MoltModel::PrintFullFemale(void)
{
  for(auto it=female.begin();it!=female.end();it++)
   {
    int female=(*it);
    string ID=to_string(female);
    string filenameFemale=pathFemale+suffix+"-animal-"+ID+".csv";
    moltfile=fopen(filenameFemale.c_str(),"a");

    set<int> herYrs=femaleyr[female];
    for(auto itf=herYrs.begin();itf!=herYrs.end();itf++)
      {
       int yr=(*itf);
       fprintf(moltfile,"%6d \t %6d \t %4d \t %6.3lf \t %6.3lf \t %6.3lf \t %6.3lf\n", 
               counter, female, yr, femalemolt[female][yr], femalemean[female], femaleslope[female][yr], femalemeanslope[female]);
      }
    fclose(moltfile);
   }
}

// Called at each step, print every year's current molt estimates.
void MoltModel::PrintFullYear(void)
{
  for(auto it=year.begin();it!=year.end();it++)
   {
    int yr=(*it);
    fprintf(yearfile,"%6d \t %4d \t %6.3lf \n", counter, yr, yearmean[yr]);
   }
}



// Following the R function, where the second parameter is the difference between x at y=.05 and y=.95. 
NumVector MoltModel::logisticMolt(NumVector d, NumVector param)
{
    double limit=0.95, middle=param[0], slope=param[1];
    double xlimit=log(limit/(1-limit));
    double b=2*xlimit/slope;
    double a=(-middle*b);
    NumVector reparam{a,b};

    Array2D<NumVector> X(d.size(),1);
    X.SetOneCol(0,d);   

    bool log_it=FALSE;   
    NumVector fit=LogisticStandard(&X,reparam,log_it);

    if(test==7 && counter%stepping[2]==0)
     {
       Test(d);
       Test(param);
       Test(reparam);
       Test(fit);
     }
    
    return(fit);
}




#endif

