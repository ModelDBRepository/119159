#ifndef _define_h_
#define _define_h_
#include<iostream.h>
#include<fstream.h>
#include<math.h>
#include <vector>
#include <algorithm>
#include<string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../../randomc.h"

#ifndef _NO_NAMESPACE
using namespace std;
#define STD std
#else
#define STD
#endif
#endif
void
Hist (vector < double >&y, vector < double >&out, int bins)
{
  for (int i = 0; i < bins; i++)
    out.push_back (0.);
  sort (y.begin (), y.end ());
  double width = y[y.size () - 1] / bins;
  for (int i = 0; i < bins; i++)
    for (unsigned int j = 0; j < y.size (); j++)
      {
	if (y[j] >= i * width && y[j] < (i + 1) * width)
	  out[i] += 1;
      }
}

void
FirstOrderStat (vector < double >&y, double &mean, double &std)
{

  double sum = 0., sumsq = 0.;
  sort (y.begin (), y.end ());
  for (unsigned int i = 0; i < y.size (); i++)
    sum += y[i];
  mean = sum / y.size ();

  for (unsigned int i = 0; i < y.size (); i++)
    sumsq += pow ((y[i] - mean), 2.0);
  std = sqrt (sumsq / (y.size () - 1));
}

void
Frate (double *t, vector < double >&y, double &r, double resolution,int filter)
{
  double sum = 0;

  for (unsigned int i = 0; i < y.size (); i++)
    {
      if (filter == 1)
	{
	  if (*t - y[i] > -resolution/2 && *t - y[i] <= resolution/2)
	    sum += 1. / resolution;
	}
      if (filter == 2)
	{
	  if (*t - y[i] > 0)
	    sum += pow((1.0/resolution),2.0)  * (*t - y[i]) * exp (-(1.0/resolution) * (*t - y[i]));
	}
      if (filter == 3)
	{
	  sum +=
	    (1. / (sqrt (2 * 3.14) * resolution)) *
	    exp (-pow ((*t - y[i]) / (resolution), 2.));
	}
    }
  r = sum;
}

void Correlate(vector<double>&x,vector<double>&y,vector <double>&C)
{
double mx,my,sigx,sigy;
double temp;

FirstOrderStat(x,mx,sigx);
FirstOrderStat(y,my,sigy);
for (unsigned int i=0;i<x.size();i++)
{
temp=0;
for (unsigned int j=0;j<x.size();j++)
{
if (i+j<y.size())
temp+=(x[j]-mx)*(y[j+i]-my)/(sigx*sigy);
}
C.push_back(temp/(x.size()-1));
}
}
