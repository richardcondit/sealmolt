
#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H 1

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
#include "/home/condit/programs/utilities/statistics.h" 

using namespace std;

// Basic likelihood functions. These are copied from likelihood.h in Rcpp folder, but with Rcpp cleaned out.

// Generic likelihood function based on Gaussian. This version accepts predicted y already calculated, and a vector of both mean and SD. Passing as pointers is huge improvement in speed. Missing values may be present, which could happen if any of data are missing. But then it is necessary to check whether every single point is missing, in which case likelihood must be -Inf.
double gaussLogLike(NumVector *predy, NumVector *y, NumVector *sd)
 {
     NumVector llike=pdfNorm(*y,*predy,*sd,1);
     int i, len=llike.size(), valid=0;
     double total=0;
     
     for(i=0;i<len;i++) if(!std::isnan(llike[i])) 
      {
          total+=llike[i];
          valid++;
       }
       
     if(valid==0) total=(-1)*INFINITY;  
    
     return(total);
 }

// Generic likelihood function based on Gaussian. This version accepts predicted y already calculated, and atomic mean and SD. 
double gaussLogLike(double mean, NumVector *y, double sd)
 {
     NumVector llike=pdfNorm(*y,mean,sd,1);
     int i, len=llike.size(), valid=0;
     double total=0;
     
     for(i=0;i<len;i++) if(!std::isnan(llike[i])) 
      {
          total+=llike[i];
          valid++;
       }
       
     if(valid==0) total=(-1)*INFINITY;  
    
     return(total);
 }

// A version accepting the observations as a map 
double gaussLogLike(double mean, map<string,double> &y, double sd)
 {
     int valid=0;
     double total=0;
     
     for(auto it=y.begin();it!=y.end();it++)  
      {
          double llike=pdfNorm(it->second,mean,sd,1);
          if(!std::isnan(llike))
           {
            total+=llike;
            valid++;
           }
       }
       
     if(valid==0) total=(-1)*INFINITY;  
     return(total);
 }

// A version accepting the observations as a map with int index
double gaussLogLike(double mean, map<int,double> &y, double sd)
 {
     int valid=0;
     double total=0;
     
     for(auto it=y.begin();it!=y.end();it++)  
      {
          double llike=pdfNorm(it->second,mean,sd,1);
          if(!std::isnan(llike))
           {
            total+=llike;
            valid++;
           }
       }
       
     if(valid==0) total=(-1)*INFINITY;  
     return(total);
 }

// Generic likelihood function based on binomial. This version accepts predicted y already calculated. Passing as pointers is huge improvement in speed. Missing values are skipped. If all values are missing, -infinity is returned. 
double binomLogLike(IntVector *y, IntVector *trials, NumVector *predprob)
 {
     NumVector llike=pdfBinom(*y,*trials,*predprob,1);
     int i, len=llike.size(), valid=0;
     double total=0;
     
     for(i=0;i<len;i++) if(!std::isnan(llike[i])) 
       {
           total+=llike[i];
           valid++;
       }
       
     if(valid==0) total=(-1)*INFINITY;  
     return(total);
 }

// Generic likelihood function based on Poisson. This version accepts predicted y already calculated. Passing as pointers is huge improvement in speed. Missing values are skipped. If all values are missing, -infinity is returned. 
double poissLogLike(IntVector *y, NumVector *pred)
 {
     NumVector llike=pdfPoisson(*y,*pred,1);
     int i, len=llike.size(), valid=0;
     double total=0;

     for(i=0;i<len;i++) 
       {
           if(std::isnan(llike[i])) return(-INFINITY);
           total+=llike[i];
           valid++;
       }
       
     if(valid==0) return(-INFINITY); 
     return(total);
 }



#endif
