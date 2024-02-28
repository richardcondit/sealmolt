/*
 A class for a 2D array of double. 
*/

#ifndef ARRAY2D_H
#define ARRAY2D_H 1

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
#include "/home/condit/programs/utilities/utilities.cpp" 
#include "/home/condit/programs/utilities/util.h"
#include "/home/condit/programs/utilities/randomgenerator.h"

 
/* A 2D array of any type T. Can be declared without arguments, or arguments to set rows and columns, or declared with R dataframe to fill, or with an existing Array2D. It can then have the size set by SetSize with rows and columns. 
 The [] operator is added so rows and individual elements can be accessed. 

Jan 2022 I need to get R and Rcpp out of this so it can be used separately. This means omitting the use of R dataframe. Instead, there will be a function in utilitiesRcpp.cpp that reads a dataframe, then uses it to fill an Array2D.
 */
 

template<typename T> class Array2D
{
  protected:
    int Ncol, Nrow;
    vector<T> A;
  public:
    Array2D();
    Array2D(int rows, int cols);
    Array2D(vector<T>);
    // Array2D(DataFrame xtable);
    
    void SetSize(int rows, int cols);
    void AddOneRow(T onerow);
    void SetOneRow(int, T);
    void SetOneCol(int j, T onecol);
    void Show();
    // void FillFromTable(vector<NumVector>);
    // void FillFromTable(DataFrame xtable);

    // Overloading the [] allows a row to be accessed from object Array2D with [], and one element with [][]; to get column requires the function 
    inline T& operator[](int i) { return((A[i])); }
    T GetOneCol(int j);
    inline IntVector Dim() { IntVector d={Nrow,Ncol}; return(d); }

    inline int rows() { return(Nrow); }
    inline int cols() { return(Ncol); }
    inline vector<T> full() { return(A); }
} ;

// Declare an empty array without setting size
template<typename T> Array2D<T>::Array2D() { }

// Or declare array with size set
template<typename T> Array2D<T>::Array2D(int rows, int cols) 
{
	SetSize(rows,cols); 
}

// Fill Array2D with a vector of vectors. By convention, each vector within the vector is a row. 
template<typename T> Array2D<T>::Array2D(vector<T> x)
{
  Ncol=x[0].size();
  Nrow=x.size();
  SetSize(Nrow,Ncol);
  for(auto i=0;i<Nrow;i++) SetOneRow(i,x[i]);  
  // Show();
}

/* I have to omit this in order for Array2D to be useful without Rcpp
Or declare passing an R dataframe to fill
template<typename T> Array2D<T>::Array2D(DataFrame xtable) 
{
  FillFromTable(xtable);
}
*/

// Setting rows and columns. This can be used after declaring empty. Or there is a declaration to create size initially.
template<typename T> void Array2D<T>::SetSize(int rows, int cols)
{
 Ncol=cols;
 Nrow=rows;
 A.resize(rows);

 typename vector<T>::iterator i;
 for(i=A.begin();i!=A.end();i++) i->resize(cols);
}

template<typename T> void Array2D<T>::SetOneRow(int row, T onerow)
{
 // for(int k=0;k<onerow.size();k++) cout << "Param " << k << " at row " << row << "=" << onerow[k] << endl;
 for(size_t j=0;j<onerow.size();j++)  A[row][j]=onerow[j];
 // for(int k=0;k<A[row].size();k++) cout << "Param " << k << " at row " << row << "=" << A[row][k] << endl;
}

/* This appears not to work */
template<typename T> void Array2D<T>::AddOneRow(T onerow) 
{ 
	Nrow++;
	A.resize(Nrow); 
	A[Nrow-1].resize(Ncol);
	SetOneRow(Nrow-1,onerow);
}

template<typename T> void Array2D<T>::SetOneCol(int j, T onecol)
{
 for(int i=0;i<(int)onecol.size();i++) A[i][j]=onecol[i];
}

template<typename T> T Array2D<T>::GetOneCol(int j)
{
 T onecol;
 for(int i=0;i<Nrow;i++) onecol.push_back(A[i][j]);
 return(onecol);
}

template<typename T> void Array2D<T>::Show() 
 { 
   for(auto i=0;i<Nrow;i++) 
    {
     for(auto j=0;j<Ncol;j++) cout << A[i][j] << ", ";
     cout << " row " << i << endl;
    }
 }

/* I have to omit this in order for Array2D to be useful without R
// Or declare array, filling it with dataframe from R
template<typename T> void Array2D<T>::FillFromTable(DataFrame xtable)
{
 int j;
 
 Ncol=xtable.size();
 
 T oneX=as<T>(xtable[0]); 
 Nrow=oneX.size();
 SetSize(Nrow,Ncol);
 
 for(j=0;j<Ncol;j++) 
  {
	 oneX=as<T>(xtable[j]);  	
	 SetOneCol(j,oneX);
  }  
 
}
*/

#endif
