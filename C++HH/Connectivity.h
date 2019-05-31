#ifndef _Connectivity_h
#define _Connectivity_h
#include <iostream.h>
#include <vector>
#include "randomc.h"
#include <time.h>
#include<vector>
#include "matrix.h"

#ifndef _NO_NAMESPACE
using namespace std;
#endif




class Connectivity{
public:
int32 seed;
int Nrow,Ncolumn;
int srow;
int bol2D,bol1D;
double P; //probability term
double **Twodarray;
double *Onedarray;
Connectivity();
~Connectivity();
Connectivity(int,int,double);
Connectivity(int,double);
void setRandomOned();
void setRandomTwod();
void setRandomTwodSym();
void setStaticOned(vector<int>&A);
void setStaticTwod(double **A);
void Print();
};
#endif
