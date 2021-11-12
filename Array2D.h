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
#include "/home/condit/programs/utilities/util.h"


/* Copied out of Rcpp/Array2D.h, but here cleaning out Rcpp. This meant first eliminating the option FillFromTable to import via R dataframe. 
 * Must substitute with input from a file.
 */

/* A 2D array of any type T. Can be declared without arguments, or arguments to set rows and columns, or declared with R dataframe to fill, or with an existing Array2D. It can then have the size set by SetSize with rows and columns. 
 * The [] operator is added so rows and individual elements can be accessed. 
 */
 
template<typename T> class Array2D
{
  protected:
    int Ncol, Nrow;
    vector<T> A;
    
  public:
    Array2D();
    Array2D(int rows, int cols);
    
    void SetSize(int rows, int cols);
    void AddOneRow(T onerow);
    void SetOneRow(int, T);
    void SetOneCol(int j, T onecol);

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
 for(size_t j=0;j<onerow.size();j++)  A[row][j]=onerow[j];
}

template<typename T> void Array2D<T>::AddOneRow(T onerow) 
{ 
	Nrow++;
	A.resize(Nrow); 
	A[Nrow-1].resize(Ncol);
	SetOneRow(Nrow-1,onerow);
}

// Assign an entire column. 
template<typename T> void Array2D<T>::SetOneCol(int j, T onecol)
{
 for(size_t i=0;i<onecol.size();i++) A[i][j]=onecol[i];
}

template<typename T> T Array2D<T>::GetOneCol(int j)
{
 T onecol;
 for(int i=0;i<Nrow;i++) onecol.push_back(A[i][j]);
 return(onecol);
}



#endif
