/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-02-27
  
--------------------------------------------------------------------------*/


#ifndef CN_MULTIFIRE_INPUTNEURON_CC
#define CN_MULTIFIRE_INPUTNEURON_CC

#include "CN_multifire_inputneuron.h"
#include "CN_inputneuron.cc"

multifire_inputneuron::multifire_inputneuron(int inlabel, tnvector<int> inpos,
					     double *inp= INPUT_p):
  inputneuron(inlabel, MF_I_IVARNO, MULTIFIRE_INPUTNEURON, inpos,
	      inp, MF_I_PNO)
{
  theE= p[2];
  t_last= -1.0;
}

multifire_inputneuron::~multifire_inputneuron()
{
}

double multifire_inputneuron::E(double *x)
{
  if (x[0] != t_last) {
    theE= p[2];
    if (fno > 0) {
      red_tx= x[0] - ((int) (x[0]/p[4]))*p[4];
      for (int i= 0; i < fno; i++)
      {
	if (((red_tx >= tb[i]) && (red_tx < te[i]))
	    || ((te[i] < tb[i]) && ((red_tx > tb[i]) || (red_tx < te[i]))))
	{
	  if ((te[i] < tb[i]) && (red_tx < te[i]))
	    red_tx+= p[4];
	  red_tx= red_tx - tb[i];
	  difft= red_tx - ((int) (red_tx/p[1]))*p[1];
	  if (difft <= p[0]) theE= p[3];
	}
      }
    }
    t_last= x[0];
  }
  return theE;
}

#endif





