#ifndef _define_h_
#define _define_h_
#include<iostream.h>
#include<iomanip.h>
#include<fstream.h>
#include<math.h>
#include <vector>
#include<string>
#include<cstring>
#include<sstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "../../Neuron.h"
#include "../../Synapse.h"
#include "../../rk4.h"
#include "../../Insynapse/Insynapse.h"
#include "../../Insynapse/Insynapse.cc"
#include "../../HHneuron/HHneuron.h"
#include "../../GABAAsynapse/GABAAsynapse.h"
#include "../../GABABsynapse/GABABsynapse.h"
#include "../../Henrysynapse/Henrysynapse.h"
#include "../../AMPAsynapse/AMPAsynapse.h"
#include "../../TwoDsynapse/TwoDsynapse.h"
#include "../../NMDAsynapse/NMDAsynapse.h"
#include "../../Ihcurrent/Ih.h"
#include "../../matrix.h"
#include "../../Calciumchannel.h"
#include "../../ITcalcium/IT.h"
#include "../../randomc.h"
#include "../../PDrule/PD.h"
#include "../../Couplingcurrent/Couplingcurrent.h"
#include <time.h>

#include <algorithm>
#include <vector>
#define N 4

#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix < double >Matrix;
#else
typedef matrix Matrix;
#endif




/*
Function Spikecount : Stores the times, in the dynamics when v>vpeak..
Basically stores the ISI for Izhikevich type models where u have resetting of the function
when it exceeds certian peak.

Function ISI, computes the average interspike interval from the data generated out of spikecount  

Function Max_vect, Min_vect computes the max and min from vector

Function Entropy : Divides the entire range of frequencies, into N bins, and computes the
probability of each bin. From that determines the std informatio entropy. The output is the
Relative entropy which is the Entropy computed divided by the max entropy=log(Nbins)

Function Complexity: Takes in the computed Relative Entropy and spits out, the Complexity
measure: Which  characterises the amount of order and disorder in the network.
*/

double MIN(double a, double b)
{
	double Temp;
	if (a<=b)
	Temp=a;
	else Temp=b;
	return Temp;

	}

double MAX(double a, double b)
{
	double Temp;
	if (a>=b)
	Temp=a;
	else Temp=b;
	return Temp;
	}

string
IntToString (int num)
{
  ostringstream myStream;	//creates an ostringstream object
  myStream << num << flush;

  /*
   * outputs the number into the string stream and then flushes
   * the buffer (makes sure the output is put into the stream)
   */

  return (myStream.str ());	//returns the string form of the stringstream object
}



string DoubleToStdStr(const double d)
{
   std::ostringstream     ostr;
   ostr << d;
   return ostr.str();
} 




double
ISI (vector < double >&y)
{
  vector < double >s;
  double average, sum;
int size;
size=y.size();  
sum = 0;
if (size>50)
{
  for (unsigned int i = size-50; i < y.size (); i++)
    s.push_back (y[i] - y[i - 1]);
  for (unsigned int i = 0; i < s.size (); i++)
    sum = sum + s[i];
  average = sum / (s.size ());
} 
else
{
cout<<"insufficient spike increase duration of time or the neuron did not spike "<<endl;
average=0;
}
 return (average);
}

double
max_vect (vector < double >&y)
{
  double val = 0;
  for (unsigned int i = 0; i < y.size (); i++)
    {
      val = max (y[i], val);

    }
  return (val);
}

double
min_vect (vector < double >&y)
{
  double val = 0;
  for (unsigned int i = 0; i < y.size (); i++)
    {
      if (val != 0)
	val = min (y[i], val);

    }
  return (val);
}



double
Spikepresent (double *tim, double *v, double threshold, double spikewidth,
	      double vpeak, int &boolean, double *spiketrigger)
{
  double spike = 0;
  if (*tim > 0 && *v > threshold && boolean == 0)
    {
      boolean = 1;
      spiketrigger[0] = *tim;
    }

  if (*tim >= spiketrigger[0] && *tim <= spiketrigger[0] + spikewidth)
    {
      spike = 1;
    }
  if (*v > vpeak)
    boolean = 0;
  return (spike);
}

void
Fread (char *y, Matrix & M, int Ntotal)
{

  ifstream fread;
  vector < double >read (Ntotal, 0);
  int row = 0;
  fread.open (y, ios::in);
  for (int i = 0; i < Ntotal; i++)
    fread >> read[i];
  if (!fread)
    {
      cerr << "error in opening" << endl;
      exit (0);
    }
  while (!fread.eof ())
    {
      for (int j = 0; j < Ntotal; j++)
	{
	  M (row, j) = read[j];

	}
      for (int j = 0; j < Ntotal; j++)
	fread >> read[j];
      row += 1;
    }

}

/****Function to compute spike times for any given spike train****/

void
spiketimes (double *tim, double *v, double Th, int &bol, vector < double >&y)
{
  if (*v > Th && bol == 0)
    y.push_back (*tim);
  if (*v > Th)
    bol = 1;
  else
    bol = 0;
}

double
Rule (double x, double a, double b, double c, double d,double strength)
{
  double out, outmax, outmax1;
  outmax = pow (a / b, a) * exp (-a);
  outmax1 = pow (c / d, c) * exp (-c);


  if (x > 0)
    out = strength * (pow (x, 1. * a) * exp (-b * fabs (x))) / outmax;
  else
    {
      if (fmod (c, 2.0) == 0)
	{
	  out = -strength * (pow (x, c) * exp (-d * fabs (x))) / (1 * outmax1);
	}
      else
	{
	  out = strength * (pow (x, c) * exp (-d * fabs (x))) / (1 * outmax1);
	}
    }
  return out;
}

double LinearRule(double x,double strength)
{
	double out;
	if (x>=-10 && x<=10)
	out=0.1*x;
	if (x>=-25 && x<-10)
	out=-(1.0/15)*x-5.0/3.0;
	if (x>10 && x<=25)
	out=-(1.0/15)*x-5.0/3.0;
	if (x>25||x<-25)
	out=0.0;
	
	out=out*strength;
	return(out);
	
	}




                  

double
box_muller (double m, double s, double &x1, double &x2, double &w)	/* normal random variate generator */
{
  /* mean m, standard deviation s */
  double y1;
  static float y2;
  static int use_last = 0;

  if (use_last)			/* use value from previous call */
    {
      y1 = y2;
      use_last = 0;
    }
  else
    {
      w = sqrt ((-2.0 * log (w)) / w);
      y1 = x1 * w;
      y2 = x2 * w;
      use_last = 1;
    }

  return (m + y1 * s);
}

double
Average (vector < double >y)
{
  double result;
  double sum = 0;
  int k = y.size ();
  if (k != 0)
    {
      for (int i = 0; i < k; i++)
	sum = sum + y[i];
      result = sum / k;
    }
  else
    result = 0.0;
  return (result);
}

void
Raster (ofstream & Outfile, double &tim, int &bul, double &volt,
	double threshold, int i)
{
//rastercomputation
  //takes in the outfile file, tim, bul: a boolean number, volt the neurons spiking, the threshold to detect spike, and spits out the time of neuron "i" firing
  bul = 0;
  if (volt > threshold && bul == 0)
    {
      Outfile << tim << " " << i  << endl;
      bul = 1;
    }
  else;
}

void
alpha_syn (double *tim, double g, double alph, vector < double >&Tspk,
	   vector < double >&sum, int &i)
{
  double *x;
  x = new double[1000];
  for (unsigned int j = 0; j < Tspk.size (); j++)
    {
      x[j] = *tim - Tspk[j];
      if (x[j] > 0)
	{
	  sum[i] = sum[i] + g * x[j] * alph * alph * exp (-alph * x[j]);
	}

    }
  delete x;
}

double
Maprule (vector < double >Tpre, vector < double >Tpost, double a, double b,
	 double c, double d)
{
  double out = 0;
  double diff;

  vector < double >Result;
  int sizeTpost = Tpost.size ();
  if (sizeTpost == 1)
    {
      for (unsigned int j = 0; j < Tpre.size (); j++)
	{
	  diff = Tpost[0] - Tpre[j];
	  Result.push_back (diff);
	}
    }

  if (sizeTpost==2||sizeTpost>2)
     {
     for (unsigned int j=0;j<Tpre.size();j++)
     {
     if (Tpre[j]<Tpost[0])
     diff=Tpost[0]-Tpre[j];
     else if(Tpre[j]>Tpost[1])
     diff=Tpost[1]-Tpre[j];
     else if (Tpre[j]>=Tpost[0] && Tpre[j]<=Tpost[1])
     {
     diff=Tpost[0]-Tpre[j];
     Result.push_back(diff);
     diff=Tpost[1]-Tpre[j]; 
     } 
     else;  
     Result.push_back(diff); 
     } 
     }
  
 /* double wold = 1000;
  if (sizeTpost >= 2)
    {
      for (int j = 1; j < Tpre.size (); j++)
	{
            wold=1000;
	  for (int i = 0; i < sizeTpost; i++)
	    {
        	      double u = Tpost[i] - Tpre[j];
	  //    cout << " " << sizeTpost << " " << Tpre.
	//	size () << " " << u << endl;
      	    
  if (fabs (u) >= 3 && fabs (u) < fabs (wold))
             wold=u;	
	else;
	    
}
	  diff = wold;
	//  cout << diff << endl;
	  Result.push_back (diff);
	}
    }
*/

  for (unsigned int j = 0; j < Result.size (); j++)
    out += Rule (Result[j], a, b, c, d,.5);
//out+=Rule(Result[j],3.0,.5,5.0,.25);

  return (out);
}

void
Pulsedetect (double *tim, vector < double >&Tspike, int &count, double &pulse)
{
  for (unsigned int j = 0; j < Tspike.size (); j++)
    {
      if (*tim > Tspike[count] + .5)
	count += 1;
      if (*tim > Tspike[count] && *tim < Tspike[count] + .5)
	pulse = 1;
      else
	pulse = 0;
    }

}

double Kuramoto(vector<double> &Tpre,vector<double> &Tpost,double period)
{
double phaseterm, costerm, sinterm, cossum = 0, sinsum = 0, kura;
double size=min(Tpre.size(),Tpost.size());
 for (int i = size-300; i < size; i++)
	{			
   phaseterm = fmod (fabs (Tpost[i] - Tpre[i]), period);
   
	  costerm = cos (phaseterm);
	  sinterm = sin (phaseterm);
	
	  cossum += costerm;
	  sinsum += sinterm;
	//cout<<phaseterm<<" "<<sinsum<<" "<<cossum<<endl;
	}
      kura = sqrt (cossum * cossum + sinsum * sinsum) / 300;
return(kura);
}
void
SlowRise (double &tim,double Threshold, double &volt, int &bul, double tauR,double &ref_time, double &out)
{
  
  if (volt > Threshold && bul==0)
     {    ref_time=tim;
      // cout<<ref_time<<endl;
        }
    if (volt>Threshold) 
    bul = 1;
else bul=0;
    if (tim >= ref_time && tim < ref_time + tauR)
    {out = 1.0;}
  else
    out = 0.0;

}
double Phase_Compute(double pre,double post,double period)
{
double phase,inter;
inter=fabs((pre-post));
phase=fmod(inter,period);
//phase=(pre-post)-floor(inter)*period;
return(phase);
}

#endif
