// Prior to Jan 2021 I had never put the ifndef line! 

#ifndef UTILITIES_CPP
#define UTILITIES_CPP


/*
Basic program utilities.  There is a long series of vector and matrix functions
for setting aside memory for a 1D vector, a 2D array (or matrix), and a 3D array.
The allocating functions have to be named differently depending on what data type
is allocated for the vector or matrix, with different versions for characters, 
integers, floats, and somestructures.  This idea was taken from Press et al.

Usage:

ivector(n) allocates memory for a vector of integers of length n;
likewise for fvector (vector of floats), cvector (char), spinfo_vector (vector
of species_information structures), and tree_vector (vector of treedata
structures). 

Likewise for imatrix, fmatrix, etc., with one called sppquad_matrix to allocate
for a matrix of structure sppquad.  Both row and column must be submitted, eg
ivector(r,c).

The deallocating functions are called free_vector and free_matrix.  Because these
are overloaded (since they take different arguments), only one name is needed for
all types.  So:

free_vector(int *v) and free_vector(float *v) both work.

There is one function for a 3D array of ints and one for floats: i3Darray(a,b,c)
and f3Darray(a,b,c).  One overloaded function, free3Darray(int ***a or float ***a)
works for both.

The paste routine is an extenstion of strcat.  By creating several with the same
name but different argument numbers, it allows up to 6 strings to be concatenated.

The strall routine searches for a single character anywhere in a string, returning
1 if it finds it, 0 if not.
*/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <pwd.h>
#include <grp.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>

// #include "/home/condit/programs/mksplist/mksplist.h"
#include "/home/condit/programs/utilities/util.h"
// #include "/home/condit/programs/diversity/diversity.h"
#include "/home/condit/programs/makefiles/plotwater.h"
// Forced to include this for the overloaded Test functions, added Nov 2021
#include "/home/condit/programs/utilities/Array2D.h"



using namespace std;

int* ivector(int n)
{
 int *x;

 x = new int[n];
 if(!x) memoryerror("could not allocate integer vector",n);

 return x;
}


void free_vector(int *v)
{
 delete[] v;
}


long* lvector(int n)
{
 long *x;

 x = new long[n];
 if(!x) memoryerror("could not allocate integer vector",n);

 return x;
}


void free_vector(long *v)
{
 delete[] v;
}


float* fvector(int n)
{
 float *x;

 x = new float[n];
 if(x==0) memoryerror("could not allocate float vector",n);

 return x;
}


void free_vector(float *v)
{
 delete[] v;
}


char* cvector(int n)
{
 char *x;

 x = new char[n];
 if(x==0) memoryerror("could not allocate character vector",n);

 return x;
}


void free_vector(char *v)
{
 delete[] v;
}

/*
struct treedata* tree_vector(int n)
{
 struct treedata *x;

 x = new treedata[n];
 if(x==0) memoryerror("could not allocate treedata structure vector",n);

 return x;
}


void free_vector(struct treedata *v)
{
 delete[] v;
}


struct species_information* spinfo_vector(int n)
{
 struct species_information *x;

 x = new species_information[n];
 if(x==0) memoryerror("could not allocate species_information structure vector",n);

 return x;
}


void free_vector(struct species_information *v)
{
 delete[] v;
}
*/

struct pointpair* pointpair_vector(int n)
{
 struct pointpair *x;

 x = new pointpair[n];
 if(x==0) memoryerror("could not allocate species_information structure vector",n);

 return x;
}


void free_vector(struct pointpair *v)
{
 delete[] v;
}


int** imatrix(int r, int c)
{
 int i, **x;

 x = new int*[r];
 for(i=0;i<r;i++) x[i] = new int[c];
 if(x==0) memoryerror("could not allocate integer matrix",r,c);

 return x;
}
   

void free_matrix(int **m, int r)
{
 int i;

 for(i=0;i<r;i++) delete[] m[i];
 delete[] m;
}


long** lmatrix(int r, int c)
{
 int i;
 long **x;

 x = new long*[r];
 for(i=0;i<r;i++) x[i] = new long[c];
 if(x==0) memoryerror("could not allocate integer matrix",r,c);

 return x;
}
   

void free_matrix(long **m, int r)
{
 int i;

 for(i=0;i<r;i++) delete[] m[i];
 delete[] m;
}


float** fmatrix(int r, int c)
{
 int i;
 float **x;

 x = new float*[r];
 for(i=0;i<r;i++) x[i] = new float[c];
 if(x==0) memoryerror("could not allocate float matrix",r,c);

 return x;
}
   

void free_matrix(float **m, int r)
{
 int i;

 for(i=0;i<r;i++) delete[] m[i];
 delete[] m;
}


char** cmatrix(int r, int c)
{
 int i;
 char **x;

 x = new char*[r];
 for(i=0;i<r;i++) x[i] = new char[c];
 if(x==0) memoryerror("could not allocate character matrix",r,c);

 return x;
}
   

void free_matrix(char **m, int r)
{
 int i;

 for(i=0;i<r;i++) delete[] m[i];
 delete[] m;
}

/*
struct sppquad** sppquad_matrix(int r, int c)
{
 int i;
 struct sppquad **x;

 x = new struct sppquad*[r];
 for(i=0;i<r;i++) x[i] = new struct sppquad[c];
 if(x==0) memoryerror("could not allocate sppquad structure matrix",r,c);

 return x;
}
   

void free_matrix(struct sppquad **m, int r)
{
 int i;

 for(i=0;i<r;i++) delete[] m[i];
 delete[] m;
}
*/

struct pointpair** pointpair_matrix(int r, int c)
{
 int i;
 struct pointpair **x;

 x = new struct pointpair*[r];
 for(i=0;i<r;i++) x[i] = new struct pointpair[c];
 if(x==0) memoryerror("could not allocate sppquad structure matrix",r,c);

 return x;
}
   

void free_matrix(struct pointpair **m, int r)
{
 int i;

 for(i=0;i<r;i++) delete[] m[i];
 delete[] m;
}


int ***i3Darray(int a, int b, int c)
{
 int i, j;
 int ***array;

 array = new int**[a];
 for(i=0;i<a;i++) array[i] = new int*[b];
 for(i=0;i<a;i++) for(j=0;j<b;j++) array[i][j] = new int[c];
 if(array==0) memoryerror("could not allocate 3D array of integers",a,b,c);

 return array;
}


void free_3Darray(int ***array, int a, int b)
{
 int i, j;

 for(i=0;i<a;i++) for(j=0;j<b;j++) delete[] array[i][j];
 for(i=0;i<a;i++) delete[] array[i];
 delete[] array;
}


float ***f3Darray(int a, int b, int c)
{
 int i, j;
 float ***array;

 array = new float**[a];
 for(i=0;i<a;i++) array[i] = new float*[b];
 for(i=0;i<a;i++) for(j=0;j<b;j++) array[i][j] = new float[c];
 if(array==0) memoryerror("could not allocate 3D array of integers",a,b,c);

 return array;
}


void free_3Darray(float ***array, int a, int b)
{
 int i, j;

 for(i=0;i<a;i++) for(j=0;j<b;j++) delete[] array[i][j];
 for(i=0;i<a;i++) delete[] array[i];
 delete[] array;
}


void memoryerror(const char *s, int n)
{
 printf("%s of size %d\n", s, n);
 exit(0);
}


void memoryerror(const char *s, int n, int m)
{
 printf("%s of size %d by %d\n", s, n, m);
 exit(0);
}


void memoryerror(const char *s, int a, int b, int c)
{
 printf("%s of size %d by %d by %d\n", s, a, b, c);
 exit(0);
}



int strall(char *s, char c)
 {
  int i, length;

  length=strlen(s);
  
  for(i=0;i<length;i++) if(s[i]==c) return 1;

  return 0;
 }


// My own substring function, returning first and last position
// of the substring s2 within string s1.

void findsubstring(char *s2, char *s1, int *position)
 {
  int i, j, len1, len2;
  char *teststr;

  position[0]=-1;
  position[1]=-1;

  len1=strlen(s1);
  len2=strlen(s2);

  teststr=cvector(len2+1);

  for(i=0;i<=(len1-len2);i++) 
   { 
    for(j=0;j<len2;j++) teststr[j]=s1[i+j];
    teststr[j]='\0';

    if(strcmp(teststr,s2)==0) 
     {
      position[0]=i;
      position[1]=i+j-1;
      break;
     }
   }

 }

// A substr function. Added Apr 2009
char *substr(char *s, int start, int end)
{
 int k, pos=0;
 char *result;
 
 result=cvector(end-start+2);
 
 for(k=start;k<=end;k++) result[pos++]=s[k];
 result[pos]='\0';
 
 return(result);
}

/*
This is the function from getline, extracted to use alone.
It takes a long string, and cuts it into shorter strings
around the char separator. Note that it always splits on tabs
as well as the submitted separator. Unlike the getstrings version,
this gives the option of treating consecutive delimiters as one
or multiple.

NEVER USED HERE AND SHOULD BE DELETED. IN THE NEW VERSION OF getline,
THE MEMBER FUNCTION IS PASSED THE SEPARATOR SO CAN BE USED INDEPENDENTLY,
OBVIATING THE NEED FOR THIS

int cutstrings(char *whole, char **pieces, char separator, bool consecutive)
  {
   int i=0, j=0, k=0, m, length;

   length=strlen(whole);
   if(whole[length-1]=='\n') m=length-2;
   else m=length-1;

   while(whole[m]==separator || whole[m]=='\t')
    {
     whole[m]='\n';
     whole[m+1]='\0';
     m--;
    }

   while(whole[i]!='\n' && whole[i]!='\0')
     {
      if(whole[i]!=separator && whole[i]!='\t') {pieces[j][k]=whole[i]; i++; k++;}
      else
       {
        if(consecutive)
         {
          if(i==0 || whole[i-1]==separator || whole[i-1]=='\t') i++;
          else {pieces[j][k]='\0'; k=0; j++; i++;}
         }
        else
         {
          if(i==0) i++;
          else {pieces[j][k]='\0'; k=0; j++; i++;}                // THIS LINE SETS HOW CONSECUTIVE SEPARATORS ARE HANDLED
         }
       }
     }
   pieces[j][k]='\0';

//   for(i=j+1;i<MAXCUTSTRING;i++) pieces[i][0]='\0';

   return(j+1);
  }
*/


// Reverses a string
char *reversestring(char *s)
{
 int i, len;
 char *bs;

 len=strlen(s);
 bs = cvector(len+1);

 for(i=0;i<len;i++) bs[i]=s[len-i-1];

 bs[len]='\0';

 return(bs);
}


// Removes lead characters (typically blanks)
char *stripleadblank(char *s, char c)
{
 int i=0, slen, start;
 char *blank;

 slen=strlen(s);

 blank=cvector(2);
 strcpy(blank,"");

 while(s[i]==c && i<slen) start=(++i);
 if(i==slen) return(blank);

 return(substr(s,start,slen-1));
}



// Removes trailing characters (typically blanks)
char *striptrailingblank(char *s, char c)
{
 int i, end;
 char *blank;

 end=strlen(s);

 blank=cvector(2);
 strcpy(blank,"");

 i=end-1;
 while(s[i]==c && i>=0) end=(--i);

 return(substr(s,0,end));
}


// Removes trailing and lead blanks from a string
char *strip(char *s, char c)
{
 return(striptrailingblank(stripleadblank(s,c),c));
}

// A version stripping lead and trailing characters, for class string
string strip(string s, string c)
{
 unsigned int start, end;
 string blank("");
 
 start=s.find_first_not_of(c);
 if(start==string::npos) return blank;
 
 end=s.find_last_not_of(c);
 
 return s.substr(start,end-start+1);
}


// Removes all instances of a char c from a string.
char *stripchar(char *s, char c)
{
 int i, k=0, slen;
 char *result;

 slen=strlen(s);
 result=cvector(slen);

 for(i=0;i<slen;i++)
  {
   if(s[i]!=c)
    {
     result[k]=s[i];
     k++;
    }
   result[k]='\0';
  }

 return(result);
}


// Removes \0 from the end of a string
char *stripbackslash0(char *s)
{
 int i, k, len;
 char *news;

 len=strlen(s);
 i=len-1;
 news=cvector(len+2);

 if(s[i]=='0' && s[i-1]=='\\')
   { for(k=0;k<len-2;k++) news[k]=s[k];  news[len-2]='\0'; }
 else strcpy(news,s);

 return(news);
}



// A tolower function for strings. Added Apr 2009
char *lowercase(char *s)
{
 unsigned int i;
 char *lower;
 
 lower=cvector(strlen(s)+2);
 
 for(i=0;i<strlen(s);i++) lower[i]=tolower(s[i]);
 lower[i]='\0';
 
 return(lower);
}

// An overloaded version of lowercase for string classes. Added Jun 2009
string lowercase(string s)
{
 string lowers=s;
 
 for(unsigned int i=0; i<s.size(); i++) lowers[i]=tolower(s[i]);
// transform(s.begin(), s.end(), lowers.begin(), ::tolower);
 return lowers;
}

// Note June 2009: with C++ strings, these paste functions are easily replaced by the + operator.
// overloaded paste functions allowing up to 7 strings to be pasted
// in one call
char *paste(char *s1, char *s2)
 {
  int len;
  char *sfinal;

  len=strlen(s1)+strlen(s2);
  sfinal=cvector(len+2);

  strcpy(sfinal,s1);
  strcat(sfinal,s2);

  return(sfinal);
 }


char *paste(char *s1, char *s2, char *s3)
 {
  return(paste(paste(s1,s2),s3));
 }


char *paste(char *s1, char *s2, char *s3, char *s4)
 {
  return(paste(paste(s1,s2,s3),s4));
 }


char *paste(char *s1, char *s2, char *s3, char *s4, char *s5)
 {
  return(paste(paste(s1,s2,s3,s4),s5));
 }


char *paste(char *s1, char *s2, char *s3, char *s4, char *s5, char *s6)
 {
  return(paste(paste(s1,s2,s3,s4,s5),s6));
 }


char *paste(char *s1, char *s2, char *s3, char *s4, char *s5, char *s6, char *s7)
 {
  return(paste(paste(s1,s2,s3,s4,s5,s6),s7));
 }

char *paste(char *s1, char *s2, char *s3, char *s4, char *s5, char *s6, char *s7, char *s8)
 {
  return(paste(paste(s1,s2,s3,s4,s5,s6,s7),s8));
 }

// From http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
void Tokenize(const string& str, vector<string>& words, const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        words.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

// My own version which can split a long string around a shorter string within
void TokenizeStr(const string origstr, vector<string>& words, const string delimiter)
{
  string str;
  str.assign(origstr);
  
  string::size_type pos = 0;
    
  while (pos != string::npos) 
   {
    pos = str.find(delimiter);
     {
      words.push_back(str.substr(0, pos));
      str.erase(0, pos + delimiter.length());
     }
   } 
  
}

// Divide a string into pieces between two delimiters. See readHTML.h, a class which invokes this on a text file. 
void TokenizeDelimPair(const string origstr, vector<string>& words, const string delim1, const string delim2)
{
  vector<string> sections;
  
  TokenizeStr(origstr, sections, delim2);
  
  if(sections.size()>0) 
   {
    vector<string>::const_iterator w;
    for(w=sections.begin(); w<sections.end()-1; w++) 
     {
      vector<string> phrases;
      string onepiece = *w;
      TokenizeStr(onepiece, phrases, delim1);
           
      if(phrases.size()>1) words.push_back(phrases.back());
     }
   }
}


// Adding functions in Dec 2020! This file may be 30 years old.
// I have some Test functions in statisticsRcpp.h, but these here so no Rcpp needed, avoiding my IntVector, NumVector.
void Test(vector<string> s) { for(auto it=s.begin();it!=s.end();it++) printf("%s ", (*it).c_str()); printf("\n"); }
void Test(vector<int> j)    { for(auto it=j.begin();it!=j.end();it++) printf("%d ", *it); printf("\n"); }
void Test(vector<long> j)    { for(auto it=j.begin();it!=j.end();it++) printf("%ld ", *it); printf("\n"); }
void Test(vector<double> x) { for(auto it=x.begin();it!=x.end();it++) printf("%lf ", *it); printf("\n"); }
void Test(Array2D<NumVector> obs) { for(auto j=0;j<obs.rows();j++) Test(obs[j]); }
void Test(Array2D<IntVector> obs) { for(auto j=0;j<obs.rows();j++) Test(obs[j]); }
void Test(map<int,double> obs) 
{
  NumVector y;
  for(auto j=obs.begin();j!=obs.end();j++) y.push_back(j->second);
  Test(y);
}
void Test(unordered_map<string,NumVector> obs) 
{
  for(auto j=obs.begin();j!=obs.end();j++) { cout << j->first << ": "; Test(j->second); }
}





#endif
