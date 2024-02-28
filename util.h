#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <map>
#include <unordered_map>


using namespace std;

#define STRLENGTH 75
#define MAXCUTSTRING 150
#define FALSE 0
#define TRUE 1

// These typedefs are used everywhere, so util.h is always needed.
typedef vector<int> IntVector;
typedef vector<double> NumVector;
typedef vector<string> StrVector;
typedef unordered_multimap<string,double> Str2Double;
typedef unordered_multimap<string,Str2Double> Str2Map2Double;
typedef unordered_multimap<string,NumVector> Str2Array;
typedef unordered_map<string,IntVector> Str2IntArray;
typedef unordered_multimap<string,int> Str2Integer;

// The rest of these are declarations for functions defined in utilities.cpp. They are ancient and not used anymore. In some of my new files, statistics.h for example, the functions are never declared. That works, providing a function is always defined before it is used.

int* ivector(int n);
long* lvector(int n);
char* cvector(int n);
float* fvector(int n);
struct treedata* tree_vector(int n);
struct species_information* spinfo_vector(int n);
struct pointpair* pointpair_vector(int n);

int** imatrix(int r, int c);
long** lmatrix(int r, int c);
float** fmatrix(int r, int c);
char** cmatrix(int r, int c);
struct sppquad** sppquad_matrix(int r, int c);
struct pointpair** pointpair_matrix(int r, int c);

int ***i3Darray(int a, int b, int c);
void free_3Darray(int ***array, int a, int b);
float ***f3Darray(int a, int b, int c);
void free_3Darray(float ***array, int a, int b);

void free_vector(int *v);
void free_vector(long *v);
void free_vector(float *v);
void free_vector(char *v);
void free_vector(struct treedata *v);
void free_vector(struct species_information *v);
void free_vector(struct pointpair *v);

void free_matrix(int **m, int r);
void free_matrix(long **m, int r);
void free_matrix(float **m, int r);
void free_matrix(char **m, int r);
void free_matrix(struct sppquad **m, int r);
void free_matrix(struct pointpair **m, int r);

void free_3Darray(int ***array, int a, int b);
void free_3Darray(float ***array, int a, int b);

void memoryerror(const char *s, int n);
void memoryerror(const char *s, int n, int m);
void memoryerror(const char *s, int a, int b, int c);

int strall(char *s, char c);
char *reversestring(char *s);
char *stripleadblank(char *s, char c);
char *striptrailingblank(char *s, char c);
char *stripchar(char *s, char c);
char *strip(char *s, char c);
string strip(string s, string c);
char *stripbackslash0(char *s);
void findsubstring(char *s2, char *s1, int *position);
char *substr(char *s, int start, int end);
char *lowercase(char *);
string lowercase(string);

char *paste(char *sfinal, char *s1);
char *paste(char *sfinal, char *s1, char *s2);
char *paste(char *sfinal, char *s1, char *s2, char *s3);
char *paste(char *sfinal, char *s1, char *s2, char *s3, char *s4);
char *paste(char *sfinal, char *s1, char *s2, char *s3, char *s4, char *s5);
char *paste(char *sfinal, char *s1, char *s2, char *s3, char *s4, char *s5, char *s6);
char *paste(char *sfinal, char *s1, char *s2, char *s3, char *s4, char *s5, char *s6, char *s7);

void Tokenize(const string& str, vector<string>& words, const string& delimiters);
void TokenizeStr(const string origstr, vector<string>& words, const string delimiter);
void TokenizeDelimPair(const string origstr, vector<string>& words, const string delim1, const string delim2);

double Sum(NumVector);

#endif

