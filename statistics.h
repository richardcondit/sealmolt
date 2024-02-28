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
#include "/home/condit/programs/utilities/util.h" 
#include "/home/condit/programs/utilities/utilities.cpp" 
#include "/home/condit/programs/utilities/Array2D.h" 

// Basic math and probability functions
// These work without preliminary declarations as long as the definition appears in the file before they are used. That is, any function that uses InvLogit has to be defined after InvLogit.

template<typename T> double Mean(map<T,double> x)
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

// Matrix multiplication, matrix %*% vector. Pass the matrix as a pointer. 
NumVector MatrixMult(Array2D<NumVector> *mat, NumVector v)
{
 NumVector result((*mat).rows(),0);
 
 if((int)v.size()!=(*mat).cols()) return(result);
 for(int j=0;j<(*mat).cols();j++) for(int i=0;i<(*mat).rows();i++) result[i]+=v[j]*(*mat)[i][j];
 
 return(result);
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

// Version 4: Gaussian probability density for a vector x, with vector for mean but scalar sd. Returns the full vector of probabilities, either logged or not. 
NumVector pdfNorm(NumVector x, NumVector mean, double sd, bool loglike)
{
    NumVector prob;
    size_t i, len=x.size();
    NumVector bad(len,NAN);
    if(mean.size()!=len) return(bad);
    
    for(i=0;i<len;i++) prob.push_back(pdfNorm(x[i],mean[i],sd,loglike));
     
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

// Negative binomial distribution for atomic observation, mean mu, and clumping parameter. Returns one probability, either logged or not. 
double pdfNegBinom(int obs, double mu, double clump, bool loglike)
{
    if(std::isnan(obs) || std::isnan(mu) || std::isnan(clump)) return(NAN);
    if(mu<=0) return(NAN);
    if(clump<=0) return(NAN);
    double p=clump/(clump+mu);
    
    NumVector oneToObs;
    for(auto i=1;i<=obs;i++) oneToObs.push_back(log((double)i));
    double logfactorial=Sum(oneToObs);

    double logProb=log(tgamma(obs+clump))-log(tgamma(clump))+clump*log(p)+obs*log(1-p)-logfactorial;
    if(!loglike) return(exp(logProb));    
    return(logProb);
}

// Negative binomial distribution for vector of observations, one mean mu, and a single clump parameter. Returns one probability, either logged or not. 
NumVector pdfNegBinom(IntVector obs, double mu, double clump, bool loglike)
{
    NumVector prob;
    size_t i, len=obs.size();
    for(i=0;i<len;i++) prob.push_back(pdfNegBinom(obs[i],mu,clump,loglike));
    return(prob);
}

// Negative binomial distribution for vector of observations, matching vector of means mu, and a single clump parameter. Returns one probability, either logged or not. 
NumVector pdfNegBinom(IntVector obs, NumVector mu, double clump, bool loglike)
{
    size_t i, len=obs.size(), lenmu=mu.size();
    NumVector prob, allNA(len,NAN);
    if(len!=lenmu) return(allNA);
    for(i=0;i<len;i++) prob.push_back(pdfNegBinom(obs[i],mu[i],clump,loglike));
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

long Min(vector<long> x)
{
    long lowest=x[0];
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

long Max(vector<long> x)
{
    long highest=x[0];
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

template<typename T> T Length(NumVector x)
{
    return((T)x.size());
}

// Calculate Gaussian kernel for y=f(x), with sd submitted
NumVector kernelGauss(NumVector x, NumVector y, double sd)
{
    int N=x.size(), i, j;
    NumVector kernel;
      for(i=0;i<N;i++)
      {
        double sumY=0, sumWt=0;
        for(j=0;j<N;j++) 
         {
           double dist=abs(x[i]-x[j]);
           double weight=pdfNorm(dist,0,sd,FALSE);
           sumY+=weight*y[j]; 
           sumWt+=weight;
         }
        kernel.push_back(sumY/sumWt);
      }
    return(kernel);
}

double InvLogit(double x)
{
 return(exp(x)/(1+exp(x)));
}

NumVector InvLogit(NumVector x)
{
 NumVector y;
 for(auto it=x.begin();it!=x.end();it++) y.push_back(InvLogit(*it));
 return(y);    
}

double Logit(double x)
{
 if(x<=0 || x>=1) return(NAN);
 return(log(x/(1-x)));    
}

NumVector Logit(NumVector x)
{
 NumVector y;
 for(auto it=x.begin();it!=x.end();it++) y.push_back(Logit(*it));
 return(y);    
}

// Repeats R logistic.standard, allowing extra parameters for asymptote and basement. Note concern over infinite or missing values, intending to be sure that missing data return prediction=nan and parameters that lead to a prediction out of range producing prediction=nan.
NumVector LogisticStandard(Array2D<NumVector> *xptr, NumVector param, bool logit)
{
 int i, Npred=(*xptr).cols(); //, N=(*xptr).rows();
 double asymp=1, basement=0, oneprob;
 NumVector X, b, pwr, y, prob, logprob;
 
 if((int)param.size()==2+Npred) asymp=param[Npred+1];
 if((int)param.size()==3+Npred) basement=param[Npred+2];    
    
 double a=param[0];
 for(i=1;i<=Npred;i++) b.push_back(param[i]);
 
 X=MatrixMult(xptr,b);
 for(auto it=X.begin();it!=X.end();it++) pwr.push_back(a+*it);
 y=InvLogit(pwr);

 for(i=0;i<(int)y.size();i++) 
    {
      if(isfinite(exp(pwr[i])) && !std::isnan(y[i]) && !std::isnan(X[i])) oneprob=y[i]*(asymp-basement)+basement;
      else oneprob=nan("");
      prob.push_back(y[i]);
      logprob.push_back(log(oneprob));
    }

 // Nov 2020 I can't get this to work. It dumpts the core if I try to use *uselog. 
 // bool *uselog = static_cast<bool*>(logptr);
 // if(*uselog) return(logprob);

 if(logit) return(logprob);
 return(prob);
}

// Repeats linear.model from R. 
NumVector LinearModel(Array2D<NumVector> *xptr, NumVector param)
{
 int i;
 NumVector y(param.size()-1,0);
 if((int)param.size()!=1+(*xptr).cols()) return(y);    

 double intercept=param[0];
 NumVector slope;
 for(i=1;i<(int)param.size();i++) slope.push_back(param[i]);
 NumVector crossprod=MatrixMult(xptr,slope);
 
 for(i=0;i<(int)crossprod.size();i++) crossprod[i]+=intercept;
 
 return(crossprod);
}


// This is linearAsymp.wt from R, but parameters are re-ordered. The break age is the third parameter; weight after senescence is fourth. The x is simply a NumVector, unlike the previous model functions
NumVector linearAsympWt(NumVector x, NumVector param)
{
 int i, N=(int)x.size();
 double senescent=log(15);
 NumVector y;
 
 double wt3=param[0], wtAsymp=param[1]+wt3, shift=param[2], wtSenescent=param[3], slope=param[1]/shift;
 for(i=0;i<N;i++) 
  {
    if(x[i]>=senescent) y.push_back(wtSenescent);
    else if(x[i]>=shift) y.push_back(wtAsymp);
    else y.push_back(wt3+x[i]*slope);  
  }
   
 return(y);
}


#endif
